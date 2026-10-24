"""Stage A3 only: repeated noisy EC, endpoint scoring, fit quality and contrasts.

Stage B3's primary comparison must use the memory epsilon rather than one-round
p_L. A future gadget plugs a RoundAdapter into memory_trial; no Stage B gadget
is built by this checkpoint. Prior approved artifacts are read-only.
"""
import argparse
from dataclasses import asdict
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

from check_order import save
from flag_design import CHECK_ORDERS,FlagDesign
from ion_noise import IdleIonNoise,IonNoise
from ion_params import H2_LIKE,JONES_MURALI_2026,DESIGN_UNIFORM_P,MEMORY_ROUNDS
from memory_fit import (fit_epsilon,compare_fits,decay_probability,
                        combine_basis_fits,R2_trend_diagnostic)
from memory_mc import (FlagMemoryRound,memory_trial,merge_memory_counts,memory_metrics,
                       memory_single_fault_checks,flag_endpoint_score)
from steane_noise import UniformNoise

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'results/stageA3_cache'
CHUNK=250_000
ORDERS=('ZZZXXX','XXXZZZ','ZXZXZX')
METRICS={'epsilon_X':'fails_X','epsilon_Z':'fails_Z'}


def one_round_product_summary(mc):
    """R=1 analogue of the headline definition, not a fit to union data."""
    n=mc['shots'];k=mc['logical_failure_counts']
    x=(k['X_only']+k['both_XZ'])/n;z=(k['Z_only']+k['both_XZ'])/n
    covariance=(k['both_XZ']/n-x*z)/n
    variance=((1-z)**2*x*(1-x)+(1-x)**2*z*(1-z))/n+2*(1-x)*(1-z)*covariance
    return {'value':1-(1-x)*(1-z),'stderr':math.sqrt(variance),
            'p_X':x,'p_Z':z,'covariance_X_Z':covariance,'definition':'R1 product summary from one-round marginals; not a union fit.'}


def fingerprint():
    names=('memory_mc.py','design_mc.py','flag_design.py','steane_ft.py','steane_noise.py',
           'ion_noise.py','ion_timing.py','ion_transport.py','ion_params.py')
    return hashlib.sha256(b'\n'.join((ROOT/'src'/p).read_bytes() for p in names)).hexdigest()


def frozen_checks():
    expected=json.loads((ROOT/'tests/fixtures/stageA2_artifacts.json').read_text())
    bad=[p for p,h in expected.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    if bad:raise RuntimeError('Approved artifact changed: '+str(bad))
    return {'artifact_count':len(expected),'all_byte_identical':True,'mismatches':[]}


def studies():
    families=[('UniformNoise(0.001)','timing-free',UniformNoise(DESIGN_UNIFORM_P[1]),1),
              ('H2_LIKE','timing-free',IonNoise(H2_LIKE),6),
              ('JONES_MURALI_2026','SEQ',IdleIonNoise(JONES_MURALI_2026),7),
              ('JONES_MURALI_2026','PAR_MR',IdleIonNoise(JONES_MURALI_2026),8)]
    for ci,(preset,policy,noise,a2ci) in enumerate(families):
        for oi,name in enumerate(ORDERS):
            design=FlagDesign(noise,CHECK_ORDERS[name],timing_policy='SEQ' if policy=='timing-free' else policy)
            yield f'f{ci}_{name}',ci,oi,preset,policy,design,f'c{a2ci}_{name}'


def fail_stop(data):
    with (ROOT/'docs/STAGE_LOG.md').open('a') as fh:
        fh.write('\n\n## Stage A3 — Memory FT/simulator failure; stopped\n\n'+json.dumps(data,indent=2)+'\n')
    raise RuntimeError('Memory fault/cross-check failure; stop for review')


def cache_point(adapter,rounds,shots,seed,path):
    config={'simulation_fingerprint':fingerprint(),'rounds':rounds,'seed':seed,
            'chunk_size':CHUNK,'numpy_version':np.__version__,'rng':'numpy.PCG64',
            'round_order':list(map(list,adapter.design.round_order)),
            'timing_policy':adapter.design.timing_policy,
            'parameters':asdict(adapter.design.params) if adapter.design.params else {'uniform_p':adapter.design.noise.p_2q}}
    if path.exists():
        data=json.loads(path.read_text())
        if data['configuration']!=config:raise ValueError(f'{path}: incompatible memory cache')
    else:data={'configuration':config,'counts':{},'rng_state':None,'chunks':0}
    rng=np.random.default_rng(seed)
    if data['rng_state'] is not None:rng.bit_generator.state=data['rng_state']
    while data['counts'].get('shots',0)<shots:
        n=min(CHUNK,shots-data['counts'].get('shots',0))
        counts=memory_trial(adapter,n,rounds,rng)
        data['counts']=merge_memory_counts(data['counts'],counts)
        data['rng_state']=rng.bit_generator.state;data['chunks']+=1
        save(path,data)
    return data


def run_mc(shots=1_000_000):
    if shots<CHUNK:raise ValueError('Use at least 250k shots per R point')
    frozen_checks();CACHE.mkdir(exist_ok=True)
    for key,ci,oi,preset,policy,d,_ in studies():
        adapter=FlagMemoryRound(d);checkpath=CACHE/f'{key}_checks.json'
        if checkpath.exists():
            checks=json.loads(checkpath.read_text())
            if checks['simulation_fingerprint']!=fingerprint():raise ValueError('Stale memory fault checks')
        else:
            checks=memory_single_fault_checks(adapter)
            # Carry every incoming single data Pauli through 16 noiseless rounds.
            xs=[];zs=[]
            for q in range(7):
                for p in 'XYZ':xs.append((1<<q) if p in 'XY' else 0);zs.append((1<<q) if p in 'ZY' else 0)
            state=(np.array(xs,np.uint16),np.array(zs,np.uint16))
            for _ in range(16):state=adapter.step(state,np.random.default_rng(0),faults={}).state
            xb,zb=flag_endpoint_score(state)
            checks.update(incoming_runs=42,incoming_failures=int(xb.sum()+zb.sum()),simulation_fingerprint=fingerprint())
            save(checkpath,checks)
        if checks['failures_union'] or checks['incoming_failures']:fail_stop({'key':key,'checks':checks})
        goal=shots
        while True:
            points=[]
            for R in MEMORY_ROUNDS:
                start=time.monotonic();path=CACHE/f'{key}_R{R}.json'
                data=cache_point(adapter,R,goal,102026+ci*10000+oi*1000+R,path)
                points.append(data['counts']);m=memory_metrics(data['counts'])
                print(key,preset,policy,'R',R,goal,'pL',m['p_L'],'seconds',round(time.monotonic()-start,2),flush=True)
            fits={label:fit_epsilon(points,metric=metric) for label,metric in METRICS.items()}
            if max(f['relative_stderr'] for f in fits.values())<=.03:break
            goal+=CHUNK;print('Extending memory points for epsilon precision',key,goal,flush=True)
        print('FITS',key,json.dumps({k:{f:x[f] for f in ('epsilon','stderr','relative_stderr','pearson_chi2','goodness_of_fit_p_value')} for k,x in fits.items()}),flush=True)


def assemble():
    frozen=frozen_checks();prior=json.loads((ROOT/'results/stageA_check_order.json').read_text())['rows']
    out={'stage':'A3','scope':'Multi-round flag memory and exact one-round pair mechanism; Stage B not started.',
         'primary_metric':'epsilon_tot','basis_fit_metrics':list(METRICS),
         'future_Stage_B3_primary_metric':'epsilon_tot=1-(1-epsilon_X)(1-epsilon_Z) from multi-round marginal fits; one-round union p_L is secondary.',
         'simulation_fingerprint':fingerprint(),
         'fit_source_sha256':hashlib.sha256((ROOT/'src/memory_fit.py').read_bytes()).hexdigest(),
         'round_lengths':list(MEMORY_ROUNDS),'fit_lengths':[2,4,8,16],
         'sampling':'Independent MC streams for every design and R; frame carried through noisy rounds; endpoint ideal decode only.',
         'union_role':'Descriptive union p_L at every R, including R=1 for continuity. The union is not fitted.',
         'trigger_attribution':'Current-round reset/measurement/2q/idle plus incoming carry, replayed as five linear components; carry is not arbitrarily reassigned to a current location type.',
         'precision_definition':'Delta-method sampling stderr from observed likelihood curvature and endpoint-binomial variances; epsilon_tot includes shared-trajectory X/Z covariance. Lack-of-fit is reported separately.',
         'approved_artifacts':frozen,'rows':{},'pairwise_comparisons':{}}
    for key,ci,oi,preset,policy,d,a2key in studies():
        checks=json.loads((CACHE/f'{key}_checks.json').read_text())
        if checks['simulation_fingerprint']!=fingerprint() or checks['failures_union'] or checks['incoming_failures']:
            raise ValueError('Missing/stale/failing memory check')
        row={'preset':preset,'policy':policy,'order':ORDERS[oi],'round_order':list(map(list,d.round_order)),
             'parameters':asdict(d.params) if d.params else {'uniform_p':d.noise.p_2q},
             'scope':'Stage 1 only (timing-free)' if preset=='H2_LIKE' else 'Noisy memory',
             'fault_checks':checks,'points':[],'fits':{},'tail_fits':{}}
        for R in MEMORY_ROUNDS:
            data=json.loads((CACHE/f'{key}_R{R}.json').read_text())
            if data['configuration']['simulation_fingerprint']!=fingerprint():raise ValueError('Stale MC point')
            m=memory_metrics(data['counts']);m['seed']=data['configuration']['seed']
            m['basis_labels']['p_L']='Descriptive logical-X-or-Z union; not fitted.'
            if d.timed:
                totals=m['round_totals'];shot_rounds=m['shots']*R
                exposure_input={'shots':shot_rounds,'triggers':totals['triggers'],
                                'no_trigger_shots':shot_rounds-totals['triggers'],
                                'trigger_circuit_counts':totals['trigger_circuit_counts']}
                exposure=d.exposure(exposure_input)
                # Round paths in one trajectory are correlated. Retain exact
                # path-weighted means, not the independent-round error formula.
                for values in exposure['monte_carlo_weighted_data'].values():
                    for metric in list(values):
                        if metric.endswith('_stderr'):values.pop(metric)
                exposure['uncertainty_note']='Means over executed rounds; independent-round exposure stderr omitted because rounds within a trajectory are correlated.'
                m['mean_idle_exposure_per_round']=exposure
            else:m['mean_idle_exposure_per_round']=None
            row['points'].append(m)
        for label,metric in METRICS.items():
            row['fits'][label]=fit_epsilon(row['points'],metric=metric)
            if row['fits'][label]['relative_stderr']>.03:raise ValueError('Epsilon needs more shots')
            row['tail_fits'][label]={str(minR):fit_epsilon(row['points'],metric=metric,minimum_rounds=minR) for minR in (4,8)}
        row['fits']['epsilon_tot']=combine_basis_fits(row['points'],row['fits']['epsilon_X'],row['fits']['epsilon_Z'])
        row['tail_fits']['epsilon_tot']={str(minR):combine_basis_fits(row['points'],row['tail_fits']['epsilon_X'][str(minR)],
                                      row['tail_fits']['epsilon_Z'][str(minR)]) for minR in (4,8)}
        row['R2_trend_diagnostics']={label:R2_trend_diagnostic(row['points'],row['tail_fits'][label]['4'],metric)
                                    for label,metric in METRICS.items()}
        row['Rge4_report_required']=any(v['outside_R2_one_sigma_MC_errorbar'] for v in row['R2_trend_diagnostics'].values())
        if row['fits']['epsilon_tot']['relative_stderr']>.03:raise ValueError('Combined metric needs more shots')
        old=prior[a2key];r1=row['points'][0]
        row['one_round_pair_attribution']={'P2':old['counts']['P2'],
             'P2_survival_corrected':old['counts']['P2_survival_corrected'],
             'P2_survival_mean_approx':old['counts']['P2_survival_mean_approx'],
             'raw_P2_over_R1_MC':old['counts']['P2']/r1['p_L'],
             'survival_P2_over_R1_MC':old['counts']['P2_survival_corrected']/r1['p_L'],
             'provenance':f'Approved A2 {a2key}, same one-round channels; this is not a multi-round P2.'}
        row['approved_A2_one_round_MC']=old['monte_carlo']
        out['rows'][key]=row
    for ci in range(4):
        comparisons=[]
        for a,b in itertools.combinations(ORDERS,2):
            ra,rb=out['rows'][f'f{ci}_{a}'],out['rows'][f'f{ci}_{b}']
            comparison={'first':a,'second':b,'direction':'second minus first','fits':{},'tail_fits':{}}
            comparison['endpoint_p_L_differences']=[]
            for pa,pb in zip(ra['points'],rb['points']):
                endpoint={'rounds':pa['rounds'],'metrics':{}}
                for label,suffix in [('union',''),('logical_X_basis_0','_X'),('logical_Z_basis_plus','_Z')]:
                    delta=pb['p_L'+suffix]-pa['p_L'+suffix]
                    error=math.hypot(pa['stderr'+suffix],pb['stderr'+suffix])
                    endpoint['metrics'][label]={'delta_p_L':delta,'combined_stderr':error,'sigma':delta/error}
                comparison['endpoint_p_L_differences'].append(endpoint)
            for metric in (*METRICS,'epsilon_tot'):
                comparison['fits'][metric]=compare_fits(rb['fits'][metric],ra['fits'][metric])
                comparison['tail_fits'][metric]={str(R):compare_fits(rb['tail_fits'][metric][str(R)],ra['tail_fits'][metric][str(R)]) for R in (4,8)}
                eps_a,eps_b=ra['fits'][metric]['epsilon'],rb['fits'][metric]['epsilon']
                sa,sb=ra['fits'][metric]['stderr'],rb['fits'][metric]['stderr']
                se=math.hypot(sb/eps_a,eps_b*sa/(eps_a*eps_a))
                rel=eps_b/eps_a-1
                comparison['fits'][metric]['relative_difference_stderr']=se
                comparison['fits'][metric]['relative_difference_ci95']=[rel-1.96*se,rel+1.96*se]
            olda,oldb=ra['approved_A2_one_round_MC'],rb['approved_A2_one_round_MC']
            d1=oldb['p_L']-olda['p_L'];se1=math.hypot(olda['stderr'],oldb['stderr'])
            comparison['approved_one_round_difference']={'delta_p_L':d1,'combined_stderr':se1,'sigma':d1/se1,
                                                        'relative_difference':d1/olda['p_L']}
            ga,gb=one_round_product_summary(olda),one_round_product_summary(oldb)
            error=math.hypot(ga['stderr'],gb['stderr']);delta=gb['value']-ga['value']
            comparison['approved_one_round_product_summary']={'first':ga,'second':gb,'delta':delta,
                  'combined_stderr':error,'sigma':delta/error,'relative_difference':gb['value']/ga['value']-1}
            comparison['approved_one_round_marginal_differences']={}
            for label,component in [('epsilon_X','X'),('epsilon_Z','Z')]:
                pa=(olda['logical_failure_counts'][component+'_only']+olda['logical_failure_counts']['both_XZ'])/olda['shots']
                pb=(oldb['logical_failure_counts'][component+'_only']+oldb['logical_failure_counts']['both_XZ'])/oldb['shots']
                err=math.sqrt(pa*(1-pa)/olda['shots']+pb*(1-pb)/oldb['shots'])
                comparison['approved_one_round_marginal_differences'][label]={'first_p':pa,'second_p':pb,
                    'delta_p':pb-pa,'combined_stderr':err,'sigma':(pb-pa)/err,'relative_difference':pb/pa-1}
            comparisons.append(comparison)
        out['pairwise_comparisons'][str(ci)]=comparisons
    mechanism_path=ROOT/'results/stageA3_pair_mechanism.json'
    if mechanism_path.exists():
        out['one_round_mechanism']={'path':'results/stageA3_pair_mechanism.json',
                                   'sha256':hashlib.sha256(mechanism_path.read_bytes()).hexdigest()}
    save(ROOT/'results/stageA3_memory.json',out);make_figures(out)
    return out


def make_figures(report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors={'ZZZXXX':'#54677d','XXXZZZ':'#c65a36','ZXZXZX':'#24858e'}
    panels=[('union','','Union endpoint (descriptive only)',None),
            ('logical_X_basis_0','_X','Logical X on |0_L>','epsilon_X'),
            ('logical_Z_basis_plus','_Z','Logical Z on |+_L>','epsilon_Z')]
    for filename,suffix,title,fitkey in panels:
        fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
        for ci,ax in enumerate(axes.flat):
            for order in ORDERS:
                row=report['rows'][f'f{ci}_{order}'];points=row['points']
                R=np.array([p['rounds'] for p in points]);y=[p['p_L'+suffix] for p in points]
                error=[1.96*p['stderr'+suffix] for p in points]
                ax.errorbar(R,y,yerr=error,fmt='o',capsize=3,color=colors[order],label=order)
                if fitkey is not None:
                    fit=row['fits'][fitkey]
                    x=np.linspace(1,16,180);ax.plot(x,decay_probability(fit['epsilon'],x),color=colors[order],alpha=.8)
            row=report['rows'][f'f{ci}_ZZZXXX']
            ax.set_title(row['preset'].replace('JONES_MURALI_2026','Jones–Murali')+' / '+row['policy'])
            ax.set(xlabel='Consecutive noisy EC rounds R',ylabel='Endpoint logical failure probability',xticks=list(MEMORY_ROUNDS))
            ax.grid(alpha=.2);ax.legend(fontsize=9)
        note='95% MC bars; no union fit' if fitkey is None else '95% MC bars; lines fit R≥2, not R=1; consult goodness of fit'
        fig.suptitle(title+'\n'+note,fontsize=13)
        fig.savefig(ROOT/f'figures/memory_rounds_{filename}.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for ci,ax in enumerate(axes.flat):
        for pos,order in enumerate(ORDERS):
            fit=report['rows'][f'f{ci}_{order}']['fits']['epsilon_tot']
            ax.errorbar(pos,fit['epsilon'],yerr=1.96*fit['stderr'],fmt='o',capsize=4,color=colors[order])
        row=report['rows'][f'f{ci}_ZZZXXX']
        ax.set_title(row['preset'].replace('JONES_MURALI_2026','Jones–Murali')+' / '+row['policy'])
        ax.set_xticks(range(3),ORDERS);ax.set_ylabel('Combined per-round epsilon_tot')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.grid(axis='y',alpha=.2);ax.set_xlim(-.4,2.4)
    fig.suptitle('Primary metric: epsilon_tot = 1 - (1-epsilon_X)(1-epsilon_Z)\n95% propagated MC bars include X/Z covariance; inspect basis fit quality',fontsize=12)
    fig.savefig(ROOT/'figures/memory_epsilon_tot.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--part',choices=('mc','report','all'),default='all')
    parser.add_argument('--shots',type=int,default=1_000_000)
    args=parser.parse_args()
    if args.part in ('mc','all'):run_mc(args.shots)
    if args.part in ('report','all'):assemble()
