"""B2-only integration validation; deliberately no performance/epsilon report."""
import json
import hashlib
import math
from pathlib import Path
import numpy as np
import stim
import steane_ft as ft
from ion_noise import IdleIonNoise
from ion_params import JONES_MURALI_2026
from check_order import save
from memory_mc import memory_trial
from prep_experiment import bernoulli
from steane_extraction import SteaneExtraction,HALF_ORDERS,syndrome,tableau_signature
from steane_memory import SteaneMemoryRound,endpoint_score
from verified_prep import run_tableau
from extraction_experiment import fingerprint,stop_failure,frozen_checks

ROOT=Path(__file__).resolve().parents[1]


def tableau_step(design,sim,faults,basis):
    log={'aborted':False,'verification':[],'syndromes':{},'corrections':{},'executed':[]}
    for hi,kind in enumerate(design.order):
        accepted=False
        for a in range(design.max_attempts):
            key=('P',hi,a);log['executed'].append(key)
            sim,bits=run_tableau(design.segments[key]['ops'],faults.get(key,{}),sim=sim)
            bit=bits['verify'];log['verification'].append({'half':kind,'attempt':a+1,'outcome':bit,'rejected':bool(bit)})
            if not bit:accepted=True;break
        if not accepted:log['aborted']=True;break
        key=('C',hi,0);log['executed'].append(key)
        sim,bits=run_tableau(design.segments[key]['ops'],faults.get(key,{}),sim=sim)
        s=syndrome(bits);fix=ft.standard_decode(s)
        log['syndromes'][kind]=s;log['corrections'][kind]=list(fix)
        for q in fix:
            if kind=='X':sim.x(q)
            else:sim.z(q)
    log['residual_signature']=tableau_signature(sim,basis)
    log['logical_failure']=bool(log['aborted'] or not ft.output_is_correctable(sim.copy(),basis))
    return sim,log


def random_faults(design,rng):
    out={}
    for key,i,op in design.locations:
        ch=design.channel((key,i,op));p=sum(ch.values())
        if rng.random()<p:out.setdefault(key,{})[i]=rng.choice(list(ch),p=np.array(list(ch.values()))/p)
    return out


def trajectory_checks(design,shots=32,rounds=4):
    adapter=SteaneMemoryRound(design);rng=np.random.default_rng(88201)
    fault_rounds=[[random_faults(design,rng) for _ in range(shots)] for r in range(rounds)]
    # Two exact, physical paths: a rejection followed by success; cap exhaustion.
    for r in range(rounds):fault_rounds[r][0]={};fault_rounds[r][1]={}
    def verifier_fault(a):
        key=('P',0,a);i=next(i for i,o in enumerate(design.segments[key]['ops']) if o[0] in ('M','MX'))
        return key,i
    key,i=verifier_fault(0);fault_rounds[0][0]={key:{i:'FLIP'}}
    fault_rounds[0][1]={key:{i:'FLIP'} for key,i in [verifier_fault(a) for a in range(4)]}
    # Leave a carried weight-one phase error after its detector in trajectory 2.
    # All faults here are explicit validation payloads, not performance samples.
    fault_rounds[0][2]={};key=('C',1,0)
    readout=next(i for i,o in enumerate(design.segments[key]['ops']) if o[0]=='IDLE' and o[1]==0 and design.segments[key]['metadata'][i]['part']=='ancilla_readout')
    fault_rounds[0][2]={key:{readout:'Z'}}
    for r in range(1,rounds):fault_rounds[r][2]={}
    encoded=[adapter.explicit_faults(f) for f in fault_rounds]
    state=adapter.initial_state(shots);frame_checks=tableau_checks=0;carried=False;rejected_seen=aborted_seen=False
    for basis in ('0','+'):
        sims=[]
        for j in range(shots):s=stim.TableauSimulator();ft.prepare_logical(s,basis);sims.append(s)
        scalar_state=[(0,0) for _ in range(shots)];dead=[False]*shots;state=adapter.initial_state(shots)
        for r in range(rounds):
            previous=state;step=adapter.step(state,rng,faults=encoded[r]);state=step.state
            if r and np.any(previous[0]|previous[1]):carried=True
            for j in range(shots):
                if dead[j]:
                    assert not state[2][j]
                    continue
                out,f=design.run(faults=fault_rounds[r][j],incoming=scalar_state[j],basis=basis)
                sims[j],t=tableau_step(design,sims[j],fault_rounds[r][j],basis)
                if f!=t or out!=(int(state[0][j]),int(state[1][j])) or bool(step.aborted[j])!=f['aborted']:
                    stop_failure({'kind':'multi_round_cross_check','order':design.half_order,'policy':design.policy,'round':r+1,'trajectory':j,'basis':basis,'frame':f,'tableau':t,'batch_frame':[int(state[0][j]),int(state[1][j])],'faults':str(fault_rounds[r][j])})
                frame_checks+=1;tableau_checks+=1;scalar_state[j]=out;dead[j]=f['aborted']
                rejected_seen |= any(v['rejected'] for v in f['verification']);aborted_seen |= f['aborted']
        counted=memory_trial(adapter,shots,rounds,rng,score=endpoint_score,fault_rounds=encoded)
        assert counted['aborted']>=1 and counted['fails_X']>=counted['aborted'] and counted['fails_Z']>=counted['aborted']
        assert counted['round_totals']['new_aborts']==counted['aborted']
    assert carried and rejected_seen and aborted_seen
    return {'trajectories_per_basis':shots,'rounds':rounds,'frame_scalar_comparisons':frame_checks,'persistent_tableau_comparisons':tableau_checks,
        'forced_rejection_seen':rejected_seen,'forced_abort_scored_in_both_bases':aborted_seen,'nonzero_frames_carried':carried,'mismatches':0,'status':'PASS'}


def run_smoke(shots=50_000):
    result={'stage':'B2','purpose':'A3 driver/Steane adapter integration smoke test only; no performance estimate, fit or comparison.',
            'simulation_fingerprint':fingerprint(),
            'harness_fingerprint':hashlib.sha256(b'\n'.join((ROOT/'src'/p).read_bytes() for p in
                ('steane_memory.py','extraction_smoke.py','memory_mc.py','design_mc.py'))).hexdigest(),
            'approved_artifacts':frozen_checks(),'rows':{},'complete':False}
    b1=json.loads((ROOT/'results/stageB1_preparation.json').read_text())
    for oi,order in enumerate(HALF_ORDERS):
        for pi,policy in enumerate(('SEQ','PAR_MR')):
            d=SteaneExtraction(IdleIonNoise(JONES_MURALI_2026),half_order=order,timing_policy=policy);adapter=SteaneMemoryRound(d)
            cross=trajectory_checks(d);key=order+'_'+policy;entries=[]
            for rounds in (1,4):
                seed=721000+oi*100+pi*10+rounds
                counts=memory_trial(adapter,shots,rounds,np.random.default_rng(seed),score=endpoint_score)
                # Each live round enters first-half first prep. Absorbed trajectories
                # stop, so this is the actual number of executed rounds, not N*R.
                executed_rounds=counts['round_totals']['segment_execution_counts'][str(('P',0,0))]
                aborts=counts['aborted'];assert aborts==counts['round_totals']['new_aborts']
                assert counts['fails_X']>=aborts and counts['fails_Z']>=aborts
                expected_time=0.;idle={str(q):0. for q in range(7)};locations={t:0 for t in ('reset','2q','measurement','idle')}
                for segkey,seg in d.segments.items():
                    n=counts['round_totals']['segment_execution_counts'][str(segkey)];res=d.resources([segkey])
                    expected_time+=n*res['duration_seconds']
                    for typ,count in res['location_counts'].items():locations[typ]+=n*count
                    for q in idle:idle[q]+=n*res['idle_exposure_per_data_qubit'][q]['seconds']
                for hi,kind in enumerate(d.order):
                    basis='plus' if kind=='X' else 'zero';ci='c7' if policy=='SEQ' else 'c8'
                    assert d.segments[('P',hi,0)]['duration']==b1['rows'][ci+'_'+basis]['duration_seconds']
                entries.append({'rounds':rounds,'trajectories':shots,'seed':seed,'status':'PASS','executed_rounds':executed_rounds,
                    'observed_aborts':bernoulli(aborts,executed_rounds),'elapsed_seconds_all_trajectories':expected_time,
                    'summed_data_idle_seconds':idle,'attempts':counts['round_totals']['attempts'],'rejections':counts['round_totals']['rejections'],
                    'executed_location_counts':locations,'segment_execution_counts':counts['round_totals']['segment_execution_counts'],
                    'rejection_cause_counts':counts['round_totals']['rejection_cause_counts'],
                    'no_crash':True,'absorbing_abort_costs_checked':True,'timing_matches_B1':True})
                print(key,'R',rounds,'PASS aborts',aborts,'/',executed_rounds,flush=True)
            result['rows'][key]={'persistent_cross_check':cross,'smoke':entries}
            save(ROOT/'results/stageB2_memory_smoke.json',result)
    result['complete']=True;save(ROOT/'results/stageB2_memory_smoke.json',result);return result


if __name__=='__main__':run_smoke()
