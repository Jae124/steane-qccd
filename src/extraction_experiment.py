"""B2 correctness checkpoint. No B3 comparisons, calibration or memory fits."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from collections import Counter
import numpy as np

import steane_ft as ft
from check_order import cases,save
from steane_extraction import SteaneExtraction,HALF_ORDERS

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'results/stageB2_cache'


def fingerprint():
    return hashlib.sha256(b'\n'.join((ROOT/'src'/p).read_bytes() for p in ('steane_extraction.py','verified_prep.py','steane_ft.py','steane_noise.py','ion_noise.py','ion_timing.py','flag_design.py','ion_params.py'))).hexdigest()


def frozen_checks():
    entries=json.loads((ROOT/'tests/fixtures/stageB1_artifacts.json').read_text())
    bad=[p for p,h in entries.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    if bad:raise RuntimeError('Approved artifacts changed: '+str(bad))
    return {'artifact_count':len(entries),'all_byte_identical':True,'mismatches':[]}


def designs():
    for ci,(preset,noise,policy) in enumerate(cases()):
        for order in HALF_ORDERS:
            yield f'c{ci}_{order}',preset,policy,SteaneExtraction(noise,half_order=order,timing_policy='SEQ' if policy=='timing-free' else policy)


def choices(design,loc,universal):
    return [v for _,v in ft.fault_list([loc[2]])] if universal else [v for v,p in design.channel(loc).items() if p>0]


def support_signature(design,universal=False):
    return hashlib.sha256(json.dumps({'order':design.order,'locations':design.locations,
        'support':[choices(design,l,universal) for l in design.locations]},sort_keys=True).encode()).hexdigest()


def stop_failure(record):
    save(ROOT/'results/stageB2_failure.json',record)
    with (ROOT/'docs/STAGE_LOG.md').open('a') as f:
        f.write('\n\n## Stage B2 — FT/cross-check failure; stopped\n\n'+json.dumps(record,indent=2)+'\n')
    raise RuntimeError('B2 failure recorded; do not tune circuit or decoder')


def check_design(design,*,universal=False):
    counts=Counter();by_part={};noiseless=[]
    def trial(faults,incoming,basis,forced,category,record,expect_abort=False):
        out,f=design.run(faults=faults,incoming=incoming,basis=basis,force_rejections=forced)
        _,t=design.run(faults=faults,incoming=incoming,basis=basis,force_rejections=forced,tableau=True)
        counts[category]+=1;counts['cross_check_runs']+=1
        if f!=t:stop_failure(dict(record,kind='frame_tableau_mismatch',basis=basis,frame=f,tableau=t,residual_frame=out))
        # Conditioned final-attempt rejection is an expected abort, not an
        # unconditional single-fault failure. All non-aborted conditioned paths
        # retain the same residual-correctability test.
        if f['aborted'] and expect_abort:counts['expected_conditioned_aborts']+=1
        elif f['logical_failure']:stop_failure(dict(record,kind='logical_failure',basis=basis,log=f,residual_frame=out))
        return f
    for basis in ('0','+'):
        f=trial({},(0,0),basis,(0,0),'noiseless_runs',{'part':'noiseless'})
        if any(f['syndromes'].values()) or f['residual_signature']!=(1,)*7:stop_failure({'kind':'noiseless','basis':basis,'log':f})
        noiseless.append({'basis':basis,'logical_preserved':True,'trivial_syndromes':True})
    for loc in design.locations:
        key,i,op=loc;meta=design.segments[key]['metadata'][i];part=('prep_'+design.segments[key]['basis']+'/'+meta['part']) if key[0]=='P' else design.order[key[1]]+'_'+meta['part']
        for payload in choices(design,loc,universal):
            record={'gadget':'SteaneExtraction','half_order':design.half_order,'policy':design.policy,'part':part,'segment':key,'op_index':i,'op':op,'payload':payload}
            for basis in ('0','+'):
                category='retry_only_fault_runs' if key[0]=='P' and key[2]>0 else 'first_attempt_or_extraction_fault_runs'
                f=trial({key:{i:payload}},(0,0),basis,(0,0),category,record)
                by_part[part]=by_part.get(part,0)+1
                if category=='retry_only_fault_runs':
                    assert key not in f['executed']
                    forced=[0,0];forced[key[1]]=key[2]
                    trial({key:{i:payload}},(0,0),basis,tuple(forced),'conditioned_retry_fault_runs',dict(record,forced_rejections=forced),expect_abort=key[2]==design.max_attempts-1)
    for q in range(7):
        for p in 'XYZ':
            incoming=((1<<q) if p in 'XY' else 0,(1<<q) if p in 'ZY' else 0)
            for basis in ('0','+'):
                trial({},incoming,basis,(0,0),'incoming_single_error_runs',{'part':'incoming','qubit':q,'payload':p})
    # Physical exhaustion: exactly four flipped verification readouts.
    # This is a multi-fault policy/cross-check test, not an FT trial.
    abort_checks=[]
    for hi in range(2):
        faults={}
        for a in range(design.max_attempts):
            key=('P',hi,a);idx=next(i for i,o in enumerate(design.segments[key]['ops']) if o[0] in ('M','MX'))
            faults[key]={idx:'FLIP'}
        for basis in ('0','+'):
            out,f=design.run(faults=faults,basis=basis);_,t=design.run(faults=faults,basis=basis,tableau=True)
            if f!=t or not f['aborted'] or not f['logical_failure']:stop_failure({'kind':'abort_policy','half':hi,'basis':basis,'frame':f,'tableau':t})
            counts['cross_check_runs']+=1;counts['physical_abort_policy_runs']+=1
            abort_checks.append({'half':hi,'basis':basis,'aborted':True,'counts_as_logical_failure':True,'resources':design.resources(f['executed'])})
    return {'catalog':'universal' if universal else 'positive preset payloads','counts':dict(counts),'fault_runs_by_part':by_part,
            'single_fault_failures':0,'incoming_failures':0,'conditioned_nonabort_failures':0,'frame_tableau_mismatches':0,
            'noiseless':noiseless,'physical_abort_checks':abort_checks,'failure_details':[],'mismatch_details':[]}


def checks():
    frozen=frozen_checks();path=ROOT/'results/stageB2_correctness.json'
    result={'stage':'B2','simulation_fingerprint':fingerprint(),'approved_artifacts':frozen,'max_attempts':4,'scope':'Full-round correctness; no B3 comparisons or epsilon fits.','unique_checks':{},'universal_checks':{},'rows':{},'complete':False}
    for key,preset,policy,d in designs():
        sig=support_signature(d)
        cache=CACHE/f'{sig}.json'
        if cache.exists():
            checked=json.loads(cache.read_text())
            if checked['simulation_fingerprint']!=fingerprint():raise ValueError('Stale B2 checks')
            checked=checked['checks']
        else:
            checked=check_design(d);save(cache,{'simulation_fingerprint':fingerprint(),'checks':checked})
        result['unique_checks'][sig]=checked
        _,log=d.run();base=d.resources(log['executed'])
        full=d.resources(list(d.segments))
        result['rows'][key]={'preset':preset,'policy':policy,'half_order':d.half_order,'half_sequence':d.order,'check_signature':sig,'checks':checked,'no_rejection_resources':base,'all_reserved_locations':full,'max_attempts':4}
        save(path,result);print(key,'PASS',checked['counts'],flush=True)
    for key,preset,policy,d in designs():
        sig=support_signature(d,universal=True)
        if sig in result['universal_checks']:continue
        cache=CACHE/f'{sig}.json'
        if cache.exists():
            cached=json.loads(cache.read_text())
            if cached['simulation_fingerprint']!=fingerprint():raise ValueError('Stale universal B2 checks')
            checked=cached['checks']
        else:
            checked=check_design(d,universal=True)
            save(cache,{'simulation_fingerprint':fingerprint(),'checks':checked})
        result['universal_checks'][sig]=checked
        save(path,result);print('Universal',d.half_order,policy,'PASS',checked['counts'],flush=True)
    result['complete']=True;save(path,result);return result


if __name__=='__main__':
    checks()
