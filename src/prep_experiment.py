"""Stage B1 checkpoint: standalone verified ancillas, no EC data coupling.

Circuit and decoder changes are forbidden after a genuine single-fault failure.
This runner saves every check before stopping. Acceptance MC is primary; exact
finite-rate XOR convolution and raw/survival-corrected P2 are benchmarks for the
standalone preparation hazard, not full-EC logical-error estimates.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from check_order import cases,save
from verified_prep import (VerifiedPreparation,BLOCK,VERIFIER,PAPER_ENCODER,PAPER_PLUS_INPUTS,
                           VERIFY_SUPPORT,CODE_GROUP,REDUCED_WEIGHT,STABILIZERS,popcount,
                           run_tableau,frame_signature,tableau_signature,weights)
import steane_ft as ft
import steane_noise as scalar

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'results/stageB1_cache'
CHUNK=1_000_000
DEFAULT_MAX_ATTEMPTS=4  # Initial preparation plus three retries.


def fingerprint():
    return hashlib.sha256(b'\n'.join((ROOT/'src'/name).read_bytes() for name in
        ('verified_prep.py','ion_noise.py','ion_timing.py','flag_design.py','steane_ft.py','steane_noise.py','ion_params.py'))).hexdigest()


def frozen_checks():
    expected=json.loads((ROOT/'tests/fixtures/stageA3_artifacts.json').read_text())
    bad=[p for p,h in expected.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    if bad:raise RuntimeError('Approved artifact changed: '+str(bad))
    return {'artifact_count':len(expected),'all_byte_identical':True,'mismatches':[]}


def gadgets():
    for ci,(preset,noise,policy) in enumerate(cases()):
        for bi,basis in enumerate(('0','+')):
            yield f'c{ci}_'+('zero' if basis=='0' else 'plus'),ci,bi,preset,policy,VerifiedPreparation(basis,noise,timing_policy='SEQ' if policy=='timing-free' else policy)


def fail_stop(row):
    with (ROOT/'docs/STAGE_LOG.md').open('a') as fh:
        fh.write('\n\n## Stage B1 — Verified-preparation FT/simulator failure; stopped\n\n'+json.dumps(row,indent=2)+'\n')
    raise RuntimeError('Verified-preparation FT/cross-check failure; stop for review')


def derivation():
    generators=[]
    for q in PAPER_PLUS_INPUTS:
        mask=1<<q
        for c,t in PAPER_ENCODER:
            if (mask>>c)&1:mask ^=1<<t
        generators.append(mask)
    group={0}
    for g in generators:group |= {v^g for v in tuple(group)}
    assert group==set(CODE_GROUP)
    verify=sum(1<<q for q in VERIFY_SUPPORT)
    assert all(popcount(verify&s)%2==0 for s in STABILIZERS) and verify not in CODE_GROUP
    g=VerifiedPreparation('0',scalar.UniformNoise(.001))
    prefix=g.base_ops[:16];cosets={};faults=0
    for idx,payload in ft.fault_list(prefix):
        if idx==7:continue # verifier reset is not an encoder fault
        x,z,_=scalar.run_ops_frame(0,0,prefix,{idx:payload});mask=(x>>7)&127
        if REDUCED_WEIGHT[mask]>1:
            representative=min((mask^s for s in CODE_GROUP),key=lambda m:(popcount(m),m))
            entry=cosets.setdefault(str(representative),{'minimal_representative':representative,'observed_masks':set(),
                                  'payload_count':0,'anticommutes_with_verifier':True})
            entry['observed_masks'].add(mask);entry['payload_count']+=1
            entry['anticommutes_with_verifier'] &= bool(popcount(mask&verify)&1)
        faults+=1
    for entry in cosets.values():entry['observed_masks']=sorted(entry['observed_masks'])
    return {'source':'https://doi.org/10.1038/srep19578','figure':'1(c)',
            'source_identifier_note':'The prompt arXiv:1610.03309 is unrelated; the correct Goto source is Sci. Rep. 6, 19578 (2016).',
            'paper_to_local_qubit_map':list(range(7)),'local_to_physical_block_map':list(BLOCK),'verifier':VERIFIER,
            'paper_encoder_CNOTs':PAPER_ENCODER,'initial_plus_local_qubits':PAPER_PLUS_INPUTS,
            'generated_X_stabilizer_supports':[[q for q in range(7) if (g>>q)&1] for g in generators],
            'same_stabilizer_row_span_as_frozen_code':True,'verification_support_local':VERIFY_SUPPORT,
            'verification_is_weight_three_logical':True,'encoder_single_fault_payloads':faults,
            'unverified_dangerous_cosets':cosets,
            'plus_derivation':'Full Hadamard conjugation: swap R/RX and M/MX, reverse every CX. No extra H primitive or unsourced timing is introduced.',
            'dangerous_component':'X for zero; Z for plus. Reduce modulo the original code stabilizers only, not the dangerous logical operator.',
            'orthogonal_component':'Canonical frame maxima are also reported modulo code stabilizers and modulo the prepared-state stabilizers, which additionally include logical Z for zero or X for plus.'}


def checks():
    frozen=frozen_checks();CACHE.mkdir(exist_ok=True)
    path=CACHE/'checks.json'
    if path.exists():
        data=json.loads(path.read_text())
        if data['simulation_fingerprint']!=fingerprint():raise ValueError('Stale B1 checks')
        if len(data.get('rows',{}))==18 and data.get('complete'):
            return data
    data={'stage':'B1','simulation_fingerprint':fingerprint(),'approved_artifacts':frozen,
          'derivation':derivation(),'universal_checks':{},'rows':{},'complete':False}
    save(path,data);seen={}
    for key,ci,bi,preset,policy,g in gadgets():
        ideal=g.noiseless_check()
        geometry=hashlib.sha256(json.dumps(g.ops).encode()).hexdigest()
        if geometry not in seen:
            universal=g.check_faults(universal=True)
            data['universal_checks'][geometry]=universal;seen[geometry]=universal
            save(path,data)
            if universal['accepted_dangerous_failures'] or universal['frame_tableau_mismatches']:
                fail_stop({'key':key,'checks':universal})
        positive=g.check_faults()
        row={'preset':preset,'policy':policy,'basis':g.basis,'noiseless':ideal,
             'geometry_signature':geometry,'checks':positive,
             'base_ops':g.base_ops,'ops':g.ops,'op_metadata':g.metadata,
             'duration_seconds':g.duration,'standalone_active_ions':8,'reserved_full_EC_layout_ions':15,
             'parameters':asdict(g.params) if g.params else {'p_2q':g.noise.p_2q},
             'scope':'Standalone 7-qubit ancilla block + verifier; data 0..6 are untouched and do not accrue B1 noise/time.'}
        data['rows'][key]=row;save(path,data)
        print(key,preset,policy,'single faults',positive['prep_fault_runs'],'accepted',positive['accepted_single_fault_runs'],
              'rejected',positive['rejected_single_fault_runs'],'failures',positive['accepted_dangerous_failures'],
              'mismatches',positive['frame_tableau_mismatches'],flush=True)
        if positive['accepted_dangerous_failures'] or positive['frame_tableau_mismatches']:fail_stop(row)
    data['complete']=True;save(path,data);return data


def merge_counts(total,addition):
    if not total:return addition
    for key,value in addition.items():
        if isinstance(value,dict):
            for label,count in value.items():total[key][label]+=count
        elif isinstance(value,list):total[key]=[a+b for a,b in zip(total[key],value)]
        elif key!='max_attempts':total[key]+=value
    return total


def sample_cache(g,shots,seed,path,*,max_attempts=None):
    config={'simulation_fingerprint':fingerprint(),'seed':seed,'chunk_size':CHUNK,'max_attempts':max_attempts,
            'numpy_version':np.__version__,'basis':g.basis,'policy':g.policy}
    if path.exists():
        data=json.loads(path.read_text())
        if data['configuration']!=config:raise ValueError('Stale B1 MC cache')
    else:data={'configuration':config,'counts':{},'rng_state':None,'chunks':0}
    rng=np.random.default_rng(seed)
    if data['rng_state'] is not None:rng.bit_generator.state=data['rng_state']
    while data['counts'].get('shots',0)<shots:
        n=min(CHUNK,shots-data['counts'].get('shots',0))
        count=(g.sample_attempts(n,rng) if max_attempts is None else g.sample_retries(n,rng,max_attempts=max_attempts))
        if max_attempts is not None:count.pop('scope');count.pop('mean_attempts')
        data['counts']=merge_counts(data['counts'],count);data['rng_state']=rng.bit_generator.state;data['chunks']+=1
        save(path,data)
    return data


def monte_carlo(shots=5_000_000,*,max_attempts=DEFAULT_MAX_ATTEMPTS,retry_caps=()):
    retry_caps=sorted(set((max_attempts,*retry_caps)))
    check=checks()
    for key,ci,bi,preset,policy,g in gadgets():
        path=CACHE/f'{key}_analysis.json'
        if path.exists():
            analysis=json.loads(path.read_text())
            if analysis['simulation_fingerprint']!=fingerprint():raise ValueError('Stale B1 analysis')
        else:
            analysis={'simulation_fingerprint':fingerprint(),'exact_finite_rate':g.exact_distribution(),'pairs':g.count_pairs()}
            save(path,analysis)
        if analysis['pairs']['P1']>0:fail_stop({'key':key,'analysis':analysis})
        data=sample_cache(g,shots,211000+ci*100+bi,CACHE/f'{key}_mc.json')
        print(key,'MC acceptance',data['counts']['accepted']/data['counts']['shots'],
              'accepted dangerous',data['counts']['accepted_dangerous'],'/',data['counts']['shots'],flush=True)
        for cap in retry_caps:
            sample_cache(g,shots,311000+ci*100+bi*10+cap,CACHE/f'{key}_attemptcap{cap}.json',max_attempts=cap)


def bernoulli(count,shots):
    p=count/shots;z=1.96
    center=(p+z*z/(2*shots))/(1+z*z/shots)
    half=z*math.sqrt(p*(1-p)/shots+z*z/(4*shots*shots))/(1+z*z/shots)
    return {'count':count,'shots':shots,'probability':p,'stderr':math.sqrt(p*(1-p)/shots),
            'ci95_Wilson':[max(0,center-half),min(1,center+half)]}


def assemble(*,max_attempts=DEFAULT_MAX_ATTEMPTS,retry_caps=()):
    retry_caps=sorted(set((max_attempts,*retry_caps)))
    frozen=frozen_checks();validated=checks()
    result={'stage':'B1','scope':'Verified preparation only; no B2 full round or B3 comparison.',
            'primary_estimate':'Monte Carlo acceptance per attempt',
            'simulation_fingerprint':fingerprint(),'approved_artifacts':frozen,'derivation':validated['derivation'],
            'universal_checks':validated['universal_checks'],'rows':{},
            'later_B3_decisions_recorded_only':{'flag_best_arm':'XXXZZZ under every preset, irrespective of prior significance',
                 'Steane_half_orders':['X-syndrome first','Z-syndrome first'],
                 'primary_metric':'epsilon_tot from A3 marginal memory fits, always report epsilon_X and epsilon_Z'},
            'max_attempts':max_attempts,
            'retry_cap_selection':'Initial attempt plus three retries: max_attempts = 4. Other caps are diagnostics only.',
            'later_B2_B3_requirements_recorded_only':{
                 'abort_accounting':'Aborts count as logical failures; report observed abort rate and its share of epsilon_tot.',
                 'cap_sensitivity':{'preset':'JONES_MURALI_2026','policy':'SEQ','max_attempts':[2,4,'infinity']},
                 'rejected_attempt_costs':'Every rejected attempt contributes real time and fault locations.'},
            'not_applicable_metrics':{'flag_trigger':'No flagged EC round in B1; report verification rejection instead.',
                                     'failure_share_through_bare':'No bare EC round in standalone preparation.'}}
    for key,ci,bi,preset,policy,g in gadgets():
        row=validated['rows'][key].copy();analysis=json.loads((CACHE/f'{key}_analysis.json').read_text())
        if analysis['simulation_fingerprint']!=fingerprint():raise ValueError('Stale B1 pair/exact analysis')
        data=json.loads((CACHE/f'{key}_mc.json').read_text());c=data['counts'];n=c['shots']
        if data['configuration']['simulation_fingerprint']!=fingerprint():raise ValueError('Stale B1 MC')
        m={'shots':n,'seed':data['configuration']['seed'],'acceptance':bernoulli(c['accepted'],n),
           'rejection':bernoulli(c['rejected'],n),'accepted_dangerous_probability':bernoulli(c['accepted_dangerous'],n),
           'dangerous_given_acceptance':bernoulli(c['accepted_dangerous'],c['accepted']),
           'verification_rejection_cause_rates':{k:v/n for k,v in c['verification_rejection_cause_counts'].items()},
           'verification_contributing_type_rates':{k:v/n for k,v in c['verification_contributing_type_counts'].items()},
           'verification_fractional_type_rates':{k:v/n for k,v in c['verification_fractional_type_counts'].items()}}
        row.update(monte_carlo=m,**{k:v for k,v in analysis.items() if k!='simulation_fingerprint'})
        row['adaptive_metrics']={'verification_trigger_rate':m['rejection']['probability'],
                                 'flag_trigger_rate':None,'failure_share_through_bare':None}
        rejection=m['rejection'];r=rejection['probability']
        row['all_attempts_fail']={
            'max_attempts':max_attempts,'per_attempt_rejection':rejection,
            'implied_from_MC_rejection':r**max_attempts,
            'propagated_stderr':max_attempts*r**(max_attempts-1)*rejection['stderr'],
            'ci95_from_rejection_Wilson':[v**max_attempts for v in rejection['ci95_Wilson']],
            'exact_finite_rate':row['exact_finite_rate']['rejection']**max_attempts,
            'scope':'One standalone preparation with independent fresh attempts; not a full-EC abort rate.'}
        hazard=m['accepted_dangerous_probability']['probability'];stderr=m['accepted_dangerous_probability']['stderr']
        row['pair_attribution_vs_MC']={k:{'estimate_over_MC':row['pairs'][k]/hazard if hazard else None,
             'estimate_minus_MC_in_stderr':(row['pairs'][k]-hazard)/stderr if stderr else None}
             for k in ('P2','P2_survival_corrected','P2_survival_mean_approx')}
        row['finite_rate_benchmark_vs_MC']={}
        for field in ('acceptance','accepted_dangerous_probability'):
            stat=m[field];diff=stat['probability']-analysis['exact_finite_rate'][field]
            row['finite_rate_benchmark_vs_MC'][field]={'MC_minus_exact':diff,'sigma':diff/stat['stderr'] if stat['stderr'] else None}
        held_out=CACHE/f'{key}_independent_validation.json'
        if held_out.exists():
            independent=json.loads(held_out.read_text())
            if independent['configuration']['simulation_fingerprint']!=fingerprint():raise ValueError('Stale held-out MC')
            counts=independent['counts'];stat=bernoulli(counts['accepted_dangerous'],counts['shots'])
            exact=row['exact_finite_rate']['accepted_dangerous_probability']
            row['independent_hazard_validation']={
                'purpose':'Fixed held-out 20M sample following low rare-hazard counts in the original H2 runs; original 5M primary rows are retained.',
                'seed':independent['configuration']['seed'],'accepted_dangerous':stat,
                'exact_probability':exact,'null_standardized_delta':(stat['probability']-exact)/math.sqrt(exact*(1-exact)/counts['shots'])}
        row['retries']={}
        for cap in retry_caps:
            cached=json.loads((CACHE/f'{key}_attemptcap{cap}.json').read_text());counts=cached['counts'];nr=counts['shots']
            if cached['configuration']['simulation_fingerprint']!=fingerprint() or cached['configuration']['max_attempts']!=cap:
                raise ValueError('Stale B1 retry cache')
            q=row['exact_finite_rate']['rejection'];a=1-q;h=row['exact_finite_rate']['accepted_dangerous_probability'];sumq=sum(q**i for i in range(cap))
            row['retries'][str(cap)]={'MC':dict(counts,mean_attempts=counts['attempt_executions']/nr,
                  exhaustion=bernoulli(counts['all_attempts_rejected'],nr)),
                'exact':{'exhaustion':q**cap,'acceptance':1-q**cap,'mean_attempts':sumq,
                         'selected_accepted_dangerous_probability':h*sumq,'standalone_failure':q**cap+h*sumq},
                'mean_duration_seconds_MC':g.duration*counts['attempt_executions']/nr if g.timed else None,
                'expected_gate_counts_MC':{'CX':11*counts['attempt_executions']/nr,'reset':8*counts['attempt_executions']/nr,'measurement':counts['attempt_executions']/nr}}
        row['idle_exposure_per_attempt']={str(q):{'seconds':math.fsum(op[2] for op in g.ops if op[0]=='IDLE' and op[1]==q),
                  'summed_Z_probability':math.fsum(-math.expm1(-op[2]/g.params.T2)/2 for op in g.ops if op[0]=='IDLE' and op[1]==q)} for q in BLOCK} if g.timed else None
        result['rows'][key]=row
    save(ROOT/'results/stageB1_preparation.json',result);draw_circuits();return result


def draw_circuits():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    fig,axes=plt.subplots(2,1,figsize=(12,8),layout='constrained')
    from verified_prep import preparation_ops
    for ax,basis in zip(axes,('0','+')):
        ops=preparation_ops(basis);ax.axvspan(7.5,15.5,color='#64748b',alpha=.09);ax.axvspan(15.5,19.5,color='#198b98',alpha=.1)
        for local,q in enumerate((*BLOCK,VERIFIER)):
            y=7-local;ax.plot([7.5,19.7],[y,y],color='#475569',lw=1)
            initial='+' if ops[local][0]=='RX' else '0';label=f'A{local} ({q})' if local<7 else f'V ({q})'
            ax.text(7.2,y,f'{label}  |{initial}>',ha='right',va='center',fontsize=10)
        for j,op in enumerate(ops[8:],8):
            if op[0]=='CX':
                ys=[7-((*BLOCK,VERIFIER).index(q)) for q in op[1:]]
                ax.plot([j,j],ys,color='#1e293b',lw=1.3);ax.plot(j,ys[0],'o',color='#1e293b',ms=5)
                ax.add_patch(Circle((j,ys[1]),.12,facecolor='white',edgecolor='#1e293b',lw=1.3));ax.plot([j-.12,j+.12],[ys[1],ys[1]],color='#1e293b',lw=1.2);ax.plot([j,j],[ys[1]-.12,ys[1]+.12],color='#1e293b',lw=1.2)
            else:ax.text(j,0,'MZ' if op[0]=='M' else 'MX',ha='center',va='center',bbox=dict(boxstyle='round,pad=.25',fc='#e9dff4',ec='#6b4696'),fontsize=10)
        ax.set_title(('|0_L> : check Z0 Z5 Z6' if basis=='0' else '|+_L> : check X0 X5 X6')+'; accept verifier bit 0',fontsize=12)
        ax.set(xlim=(5.2,20),ylim=(-.55,7.65),yticks=[],xticks=range(8,20),xlabel='Original operation index (0–7 are the eight initial resets)')
        ax.spines[['top','right','left']].set_visible(False)
    fig.suptitle('Standalone verified ancilla preparation: Goto Fig. 1(c) and Hadamard dual\nLocal A0…A6 are physical 7…13; V is 14. Circuit order, not a hardware schedule.',fontsize=13)
    fig.savefig(ROOT/'figures/verified_prep_circuits.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--part',choices=('checks','mc','report','all'),default='all')
    parser.add_argument('--shots',type=int,default=5_000_000)
    parser.add_argument('--max-attempts',type=int,default=DEFAULT_MAX_ATTEMPTS,
                        help='Total attempts, including the initial attempt (default: 4).')
    parser.add_argument('--retry-caps',type=int,nargs='*',default=[],help='Additional diagnostic total-attempt caps.')
    args=parser.parse_args()
    if args.max_attempts<1 or any(k<1 for k in args.retry_caps):parser.error('Attempt caps must be positive integers.')
    if args.part in ('checks','all'):checks()
    if args.part in ('mc','all'):monte_carlo(args.shots,max_attempts=args.max_attempts,retry_caps=args.retry_caps)
    if args.part in ('report','all'):assemble(max_attempts=args.max_attempts,retry_caps=args.retry_caps)
