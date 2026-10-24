import json
import math
import numpy as np
import pytest
import stim

import steane_ft as ft
import steane_noise as scalar
from check_order import cases, geometry_signature
from design_mc import BatchSampler, CATEGORIES, merge_counts, metrics
from flag_design import CHECK_ORDERS, FlagDesign, ec_round, ec_round_frame
from ion_noise import IdleIonNoise, _state_signature, _frame_signature
from ion_params import JONES_MURALI_2026


def scalar_oracle(sampler, shots, faults, *, tableau_limit=0):
    out=sampler.simulate(shots,faults,traces=True);t=out.pop('traces')
    per_shot={}
    for li,(indices,codes) in faults.items():
        loc=sampler.locations[li];payloads=list(sampler.channels[li])
        for shot,code in zip(indices,codes):
            per_shot.setdefault(shot,{}).setdefault(loc[:2],{})[loc[2]]=payloads[code-1]
    for j,shot in enumerate(t['shot_indices']):
        log={};d=sampler.design
        x,z=ec_round_frame(0,0,per_shot[shot],circuits=d.circuits,round_order=d.round_order,log=log)
        assert (x,z)==(t['fx'][j],t['fz'][j])
        assert scalar.logical_failure(x,z)==(t['x_fail'][j],t['z_fail'][j])
        assert (log['trigger_index'] if log['trigger_index'] is not None else -1)==t['stop'][j]
        if log['trigger']:
            assert log['trigger_bits']=={'flag':t['flag'][j],'synd':t['trigger_synd'][j]}
            assert (log['s_x'],log['s_z'])==(t['s_x'][j],t['s_z'][j])
            assert sum(1<<q for q in log['x_fix'])==t['x_fix'][j]
            assert sum(1<<q for q in log['z_fix'])==t['z_fix'][j]
            assert set(scalar.trigger_contributors(per_shot[shot],d.circuits,log))=={
                c for i,c in enumerate(CATEGORIES) if t['contributors'][i,j]}
        if j<tableau_limit:
            for basis in ('0','+'):
                sim=stim.TableauSimulator();ft.prepare_logical(sim,basis)
                assert ec_round(sim,faults=per_shot[shot],circuits=d.circuits,round_order=d.round_order)==log
                assert _state_signature(sim,basis)==_frame_signature(x,z,basis)
    assert sum(out['trigger_cause_counts'].values())==out['triggers']
    assert sum(out['logical_failure_counts'].values())==out['fails']
    assert math.isclose(sum(out['trigger_fractional_type_counts'].values()),out['triggers'])
    return out


@pytest.mark.parametrize('order',CHECK_ORDERS)
def test_batch_exhaustive_payloads_and_triggered_bare_against_scalar(order):
    seen=set()
    for _,noise,policy in cases():
        d=FlagDesign(noise,CHECK_ORDERS[order],timing_policy='SEQ' if policy=='timing-free' else policy)
        signature=geometry_signature(d)
        if signature in seen:continue
        seen.add(signature);b=BatchSampler(d);faults={};shot=0
        trigger=next(i for i,l in enumerate(b.locations) if l[:2]==('F',0) and l[3][0]=='M')
        flip=list(b.channels[trigger]).index('FLIP')+1
        for li,ch in enumerate(b.channels):
            for code,prob in enumerate(ch.values(),1):
                if not prob:continue
                faults.setdefault(li,([],[]))[0].append(shot);faults[li][1].append(code);shot+=1
                if b.locations[li][0]=='B':
                    for loc,val in ((li,code),(trigger,flip)):
                        faults.setdefault(loc,([],[]))[0].append(shot);faults[loc][1].append(val)
                    shot+=1
        faults={li:(np.array(a),np.array(c)) for li,(a,c) in faults.items()}
        scalar_oracle(b,shot+1,faults)


@pytest.mark.parametrize('policy',['SEQ','PAR_MR'])
@pytest.mark.parametrize('order',CHECK_ORDERS)
def test_batch_mixed_random_faults_and_tableau(policy,order):
    b=BatchSampler(FlagDesign(IdleIonNoise(JONES_MURALI_2026),CHECK_ORDERS[order],timing_policy=policy))
    rng=np.random.default_rng(4826)
    # High-occupancy explicit injections exercise cancellations and mixed causes,
    # independent of the production sampler and without changing gadget channels.
    faults={}
    for li,ch in enumerate(b.channels):
        positive=np.flatnonzero(np.array(list(ch.values()))>0)+1
        if len(positive):
            hits=np.flatnonzero(rng.random(600)<0.035)
            faults[li]=(hits,rng.choice(positive,len(hits)))
    scalar_oracle(b,600,faults,tableau_limit=50)


def test_sparse_sampling_location_independence_payload_weights_and_resume():
    b=BatchSampler(FlagDesign(scalar.UniformNoise(0.1)))
    n=100000;rng=np.random.default_rng(8226);faults=b.sample_faults(n,rng)
    for li,(hits,codes) in faults.items():
        assert len(np.unique(hits))==len(hits)
        assert abs(len(hits)-n*0.1)<7*math.sqrt(n*0.1*0.9)
        if b.locations[li][3][0]=='CX':
            counts=np.bincount(codes,minlength=16)[1:]
            assert np.max(np.abs(counts-len(hits)/15))<7*math.sqrt(len(hits)/15)
    a,c=[set(faults[i][0]) for i in (0,1)]
    assert abs(len(a&c)-n*0.01)<7*math.sqrt(n*0.01*0.99)
    state=json.loads(json.dumps(rng.bit_generator.state));first=b.sample(10000,rng)
    rng2=np.random.default_rng();rng2.bit_generator.state=state
    assert first==b.sample(10000,rng2)
    merged=merge_counts({},first);merged=merge_counts(merged,b.sample(10000,rng2))
    assert metrics(merged)['shots']==20000
    empty=b.simulate(100,{})
    assert metrics(empty)['p_L']==0 and empty['triggers']==0


def test_batch_rejects_synthetic_transport():
    noise=IdleIonNoise(JONES_MURALI_2026,transport_events=[('F',0,0,0,'shuttle')])
    with pytest.raises(ValueError,match='empty transport'):BatchSampler(FlagDesign(noise))


def test_sparse_sampling_unequal_payload_probabilities():
    locations=scalar.all_locations()
    selected=next(l for l in locations if l[3][0]=='CX')
    class Unequal:
        def channel(self,loc):return {'XI':0.1,'IZ':0.2} if loc==selected else {}
    b=BatchSampler(FlagDesign(Unequal()))
    fs=b.sample_faults(100000,np.random.default_rng(9926))
    assert len(fs)==1
    _,codes=next(iter(fs.values()))
    for code,p in [(1,0.1),(2,0.2)]:
        assert abs(np.count_nonzero(codes==code)-100000*p)<7*math.sqrt(100000*p*(1-p))
