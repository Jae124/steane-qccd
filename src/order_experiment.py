"""Stage A2 only: cached order MC, exact pair attribution, and MC calibration.

Every chunk records accumulated integer counts and the NumPy RNG state. Reuse
requires matching simulation/configuration fingerprints. Old artifacts are read
only. No Steane-style extraction is implemented here.
"""
import argparse
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from check_order import cases, frozen_artifacts, geometry_signature, save, source_fingerprint
from design_mc import BatchSampler, merge_counts, metrics
from flag_design import CHECK_ORDERS, FlagDesign
from ion_analysis import fingerprint as approved_fingerprint
from ion_calibration import calibrate_targets
from ion_params import DESIGN_CALIBRATION_EXTENSION_P
from steane_noise import UniformNoise, all_locations, validated_channel

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'results/stageA2_cache'
CHUNK=1_000_000


def simulation_fingerprint():
    names=('design_mc.py','flag_design.py','steane_ft.py','steane_noise.py',
           'ion_noise.py','ion_timing.py','ion_transport.py','ion_params.py')
    return hashlib.sha256(b'\n'.join((ROOT/'src'/n).read_bytes() for n in names)).hexdigest()


def studies():
    for case_index,(preset,noise,policy) in enumerate(cases()):
        for order_index,(name,order) in enumerate(CHECK_ORDERS.items()):
            design=FlagDesign(noise,order,timing_policy='SEQ' if policy=='timing-free' else policy)
            yield f'c{case_index}_{name}',case_index,order_index,preset,policy,design


def frozen_a1():
    expected=json.loads((ROOT/'tests/fixtures/stageA1_artifacts.json').read_text())
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in expected.items()}


def guard_frozen():
    if not frozen_artifacts()['all_byte_identical'] or not all(frozen_a1().values()):
        raise RuntimeError('Approved artifacts changed; stop')


def load_cache(path,configuration,default):
    if not path.exists():return dict(configuration=configuration,**default)
    data=json.loads(path.read_text())
    if data['configuration']!=configuration:
        raise ValueError(f'{path}: incompatible cache; choose a new namespace')
    return data


def fail_stop(row):
    with (ROOT/'docs/STAGE_LOG.md').open('a') as fh:
        fh.write('\n\n## Stage A2 — Fault/simulator check failure; stopped\n\n'+json.dumps(row,indent=2)+'\n')
    raise RuntimeError('Fault/simulator check failure; stop for review')


def _scale_pair_counts(base,p,base_p,active_locations):
    """Exact for equal-p active locations, including zero-rate spectators."""
    out=deepcopy(base);scale=(p/base_p)**2
    survival=(1-p)**(active_locations-2)
    for key in ('P2','P2_by_logical_type','P2_by_stage_pair','P2_by_type_pair'):
        value=base[key]
        out[key]={k:v*scale for k,v in value.items()} if isinstance(value,dict) else value*scale
    out['P2_survival_corrected']=out['P2']*survival
    out['P2_survival_by_stage_pair']={k:v*survival for k,v in out['P2_by_stage_pair'].items()}
    out['P2_survival_by_type_pair']={k:v*survival for k,v in out['P2_by_type_pair'].items()}
    out['mean_location_fault_probability']=p*active_locations/out['n_locations']
    out['P2_survival_mean_approx']=out['P2']*(1-out['mean_location_fault_probability'])**(out['n_locations']-2)
    return out


def pair_attribution():
    guard_frozen();CACHE.mkdir(exist_ok=True)
    path=CACHE/'pairs.json'
    config={'simulation_fingerprint':simulation_fingerprint(),'stage':'A2','empty_transport':True}
    data=load_cache(path,config,{'rows':{}})
    matrix=json.loads((ROOT/'results/stageA1_fault_checks.json').read_text())
    if matrix['source_fingerprint']!=source_fingerprint():raise ValueError('Stale approved A1 correctness matrix')
    reuse={}
    for key,ci,oi,preset,policy,d in studies():
        geom=geometry_signature(d)
        if key in data['rows']:
            row=data['rows'][key]
        else:
            matching=[r for r in matrix['rows'] if r['geometry_signature']==geom and r['preset']==preset]
            if len(matching)!=1:raise ValueError('Missing approved A1 fault check')
            checks=matching[0]['checks']
            if checks['single_fault_failures'] or checks['mismatches']:fail_stop(matching[0])
            if ci in (7,8) and oi==0:
                ref=json.loads((ROOT/f'results/stageA1_reference_{policy}.json').read_text())
                counts=ref['counts'];provenance=f'Approved A1 {policy} baseline counts, exact same circuits/channels.'
            elif geom in reuse:
                base,base_p,active=reuse[geom]
                counts=_scale_pair_counts(base,d.noise.p_2q,base_p,active)
                provenance='Exact p^2 reweighting of identical equal-p active payloads; survival product recomputed.'
            else:
                start=time.monotonic();print('Pair enumeration:',key,preset,policy,flush=True)
                counts=d.count_faults();provenance='Exact exhaustive pair enumeration.'
                print('  P2',counts['P2'],'seconds',round(time.monotonic()-start,2),flush=True)
            if counts['single_fault_failures']:fail_stop({'key':key,'counts':counts})
            row={'order_name':list(CHECK_ORDERS)[oi],'round_order':d.round_order,'preset':preset,'timing_policy':policy,
                 'scope':'Stage 1 only' if ci==6 else 'One adaptive EC round',
                 'counts':counts,'fault_checks':checks,'pair_provenance':provenance,
                 'p_2q':d.params.p_2q if d.params else d.noise.p_2q,
                 'parameters':asdict(d.params) if d.params else {'p_2q':d.noise.p_2q},
                 'geometry_signature':geom}
            data['rows'][key]=row;save(path,data)
        if ci<6:
            active=sum(sum(validated_channel(d,l).values())>0 for l in all_locations(d.circuits))
            reuse.setdefault(geom,(row['counts'],d.noise.p_2q,active))
    return data


def sample_cached(path,design,shots,seed):
    config={'simulation_fingerprint':simulation_fingerprint(),'seed':seed,'chunk_size':CHUNK,
            'rng':'numpy.PCG64','numpy_version':np.__version__,
            'order':list(map(list,design.round_order)),'timing_policy':design.timing_policy,
            'ordered_channel_hash':hashlib.sha256(json.dumps(
                [[l,list(validated_channel(design,l).items())] for l in all_locations(design.circuits)]).encode()).hexdigest()}
    data=load_cache(path,config,{'counts':{},'chunks':0,'rng_state':None})
    if data['counts'].get('shots',0)>=shots:return metrics(data['counts'])
    rng=np.random.default_rng(seed)
    if data['rng_state'] is not None:rng.bit_generator.state=data['rng_state']
    sampler=BatchSampler(design)
    while data['counts'].get('shots',0)<shots:
        count=min(CHUNK,shots-data['counts'].get('shots',0))
        merge= sampler.sample(count,rng)
        data['counts']=merge_counts(data['counts'],merge)
        data['chunks']+=1;data['rng_state']=rng.bit_generator.state
        save(path,data)
    return metrics(data['counts'])


def monte_carlo():
    pairs=pair_attribution();goals={}
    for ci in range(9):
        rates=[r['counts']['P2'] for key,r in pairs['rows'].items() if key.startswith(f'c{ci}_')]
        # Fixed allocation from deterministic P2, with a factor-two variance margin.
        # MC determines the estimates and the final precision check below.
        goals[ci]=max(5_000_000,CHUNK*math.ceil(4/(0.02**2*min(rates))/CHUNK))
    by_case={}
    for key,ci,oi,preset,policy,d in studies():by_case.setdefault(ci,[]).append((key,oi,preset,policy,d))
    for ci,group in by_case.items():
        goal=goals[ci]
        while True:
            rows={}
            for key,oi,preset,policy,d in group:
                start=time.monotonic()
                rows[key]=sample_cached(CACHE/f'{key}_mc.json',d,goal,83000+ci*100+oi)
                print(key,preset,policy,rows[key]['shots'],'shots',rows[key]['p_L'],'+/-',rows[key]['stderr'],
                      'seconds',round(time.monotonic()-start,2),flush=True)
            base=rows[group[0][0]]
            relative=[math.hypot(base['stderr'],r['stderr'])/min(base['p_L'],r['p_L']) for r in rows.values()]
            if max(relative)<=0.02:break
            goal+=CHUNK
            print('Extending all arms to satisfy combined MC precision:',goal,flush=True)
    return pairs


def calibration():
    guard_frozen();CACHE.mkdir(exist_ok=True)
    base=json.loads((ROOT/'results/uniform_mc_calibration.json').read_text())
    if base['configuration']['simulation_fingerprint']!=approved_fingerprint():raise ValueError('Stale approved uniform calibration')
    path=CACHE/'calibration_extension.json'
    data=load_cache(path,{'simulation_fingerprint':simulation_fingerprint(),
                          'grid':list(DESIGN_CALIBRATION_EXTENSION_P),'shots_per_point':5_000_000}, {'points':{}})
    A=json.loads((ROOT/'results/stage0_channels.json').read_text())['A']['any']
    for i,p in enumerate(DESIGN_CALIBRATION_EXTENSION_P):
        if str(p) in data['points']:continue
        design=FlagDesign(UniformNoise(p))
        mc=sample_cached(CACHE/f'grid_{p:g}.json',design,5_000_000,94000+i)
        data['points'][str(p)]={'p':p,'monte_carlo':mc,'P2':A*p*p,
                                'P2_survival_corrected':A*p*p*(1-p)**94}
        save(path,data);print('Calibration extension',p,mc['p_L'],flush=True)
    return list(base['points'].values())+list(data['points'].values())


def report():
    guard_frozen()
    pairs=pair_attribution();grid=calibration()
    result={'stage':'A2','primary_estimate':'Monte Carlo','simulation_fingerprint':simulation_fingerprint(),
            'MC_method':'Exact sparse independent location sampling; binomial count then uniform subset; NumPy PCG64; 1M-shot resumable chunks.',
            'shot_allocation':'Per family max(5M, ceil_to_1M(4/(0.02^2*min_order_P2))); extend complete family if measured combined stderr exceeds 2% of smaller p_L.',
            'difference_convention':'Order minus ZZZXXX under the same preset and timing policy; independent target MC streams.',
            'significance_rule':'abs(delta)/combined_stderr >= 3; per-comparison, no family-wise multiplicity claim.',
            'frozen_artifacts':frozen_artifacts(),'approved_A1_artifacts_unchanged':frozen_a1(),
            'uniform_calibration':{'approved_source':'results/uniform_mc_calibration.json',
                                   'approved_source_sha256':hashlib.sha256((ROOT/'results/uniform_mc_calibration.json').read_bytes()).hexdigest(),
                                   'extension_source':'results/stageA2_cache/calibration_extension.json','points':sorted(grid,key=lambda r:r['p'])},
            'rows':{}}
    targets={}
    for key,ci,oi,preset,policy,d in studies():
        data=json.loads((CACHE/f'{key}_mc.json').read_text())
        if data['configuration']['simulation_fingerprint']!=simulation_fingerprint():raise ValueError('Stale MC cache')
        row=deepcopy(pairs['rows'][key]);mc=metrics(data['counts']);row['monte_carlo']=mc
        row['MC_seed']=data['configuration']['seed'];row['MC_cache']=f'results/stageA2_cache/{key}_mc.json'
        row['logical_failure_inclusive_rates']={
            label:(mc['logical_failure_counts'][label+'_only']+mc['logical_failure_counts']['both_XZ'])/mc['shots']
            for label in ('X','Z')}
        row['idle_exposure']=d.exposure(mc)
        row['trigger_position_labels']=[f'{i}:{kind}{k}' for i,(kind,k) in enumerate(d.round_order)]
        row['trigger_position_rates']=[n/mc['shots'] for n in mc['trigger_circuit_counts']]
        targets[key]={'monte_carlo':mc,'p_2q':row['p_2q']}
        result['rows'][key]=row
    cal=calibrate_targets(grid,targets,bootstrap_draws=20000,seed=503026)
    result['calibration_statistics']=cal
    result['uniform_pseudothreshold']=cal['pseudothreshold']
    A=json.loads((ROOT/'results/stage0_channels.json').read_text())['A']['any']
    for key,row in result['rows'].items():
        ci=int(key.split('_')[0][1:]);base=result['rows'][f'c{ci}_ZZZXXX']
        mc=row['monte_carlo'];bmc=base['monte_carlo']
        error=0.0 if key.endswith('_ZZZXXX') else math.hypot(mc['stderr'],bmc['stderr'])
        delta=mc['p_L']-bmc['p_L'];relative=error/min(mc['p_L'],bmc['p_L'])
        if relative>0.02:raise ValueError(f'{key}: insufficient combined MC precision')
        row['versus_same_policy_baseline']={'delta_p_L':delta,'combined_stderr':error,
            'sigma':delta/error if error else 0,'significant_at_3sigma':abs(delta)>=3*error if error else False,
            'relative_combined_stderr':relative,'delta_P2':row['counts']['P2']-base['counts']['P2'],
            'delta_P2_survival_corrected':row['counts']['P2_survival_corrected']-base['counts']['P2_survival_corrected']}
        row['p_eff']=cal['targets'][key]
        row['p_eff_P2_secondary']=math.sqrt(row['counts']['P2']/A)
        row['pair_attribution_vs_MC']={k:{'estimate_over_MC':row['counts'][k]/mc['p_L'],
            'estimate_minus_MC_in_stderr':(row['counts'][k]-mc['p_L'])/mc['stderr']}
            for k in ('P2','P2_survival_corrected','P2_survival_mean_approx')}
    save(ROOT/'results/stageA_check_order.json',result)
    make_figure(result)
    return result


def make_figure(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names=list(CHECK_ORDERS);colors=['#66717e','#c45536','#287f8e','#8e5da2']
    fig,axes=plt.subplots(3,3,figsize=(13,11),layout='constrained')
    for ci,ax in enumerate(axes.flat):
        rows=[result['rows'][f'c{ci}_{name}'] for name in names]
        ax.errorbar(range(4),[r['monte_carlo']['p_L'] for r in rows],
                    yerr=[1.96*r['monte_carlo']['stderr'] for r in rows],fmt='none',ecolor='#334155',capsize=4)
        ax.scatter(range(4),[r['monte_carlo']['p_L'] for r in rows],c=colors,s=42,zorder=3)
        title=rows[0]['preset'].replace('JONES_MURALI_2026','Jones–Murali')
        title+='\n'+('Stage 1 only; timing-free' if ci==6 else rows[0]['timing_policy'])
        ax.set_title(title,fontsize=11)
        ax.set_xticks(range(4),names,rotation=20);ax.set_ylabel('Logical failure probability p_L')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.grid(axis='y',alpha=.25)
        ax.set_xlim(-.35,3.35)
    fig.suptitle('Stabilizer-check order: Monte Carlo estimates with 95% sampling bars\nEach panel uses one preset and one timing policy',fontsize=14)
    fig.savefig(ROOT/'figures/check_order_comparison.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--part',choices=('pairs','mc','calibration','report','all'),default='all')
    args=parser.parse_args()
    if args.part in ('pairs',):pair_attribution()
    if args.part in ('mc','all'):monte_carlo()
    if args.part in ('calibration','all'):calibration()
    if args.part in ('report','all'):report()
