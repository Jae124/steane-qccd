import math
import numpy as np
import pytest
import steane_noise as scalar
from ion_noise import IdleIonNoise
from ion_params import JONES_MURALI_2026,H2_LIKE
from ion_noise import IonNoise
from steane_extraction import (SteaneExtraction,HALF_ORDERS,DATA,BLOCK,syndrome,frame_signature)
from extraction_experiment import frozen_checks


@pytest.mark.parametrize('order',HALF_ORDERS)
@pytest.mark.parametrize('policy',['SEQ','PAR_MR'])
def test_full_round_noiseless_and_every_incoming_single_error(order,policy):
    d=SteaneExtraction(IdleIonNoise(JONES_MURALI_2026),half_order=order,timing_policy=policy)
    for basis in ('0','+'):
        for x,z in [(0,0)]+[(1<<q if p in 'XY' else 0,1<<q if p in 'ZY' else 0) for q in DATA for p in 'XYZ']:
            _,f=d.run(basis=basis,incoming=(x,z));_,t=d.run(basis=basis,incoming=(x,z),tableau=True)
            assert f==t and not f['logical_failure']
            assert f['residual_signature']==(1,)*7


@pytest.mark.parametrize('order',HALF_ORDERS)
@pytest.mark.parametrize('policy',['SEQ','PAR_MR'])
def test_all_attempts_charge_data_wait_and_transversal_live_windows(order,policy):
    d=SteaneExtraction(IdleIonNoise(JONES_MURALI_2026),half_order=order,timing_policy=policy)
    _,base=d.run();_,retried=d.run(force_rejections=(3,3))
    r0=d.resources(base['executed']);r4=d.resources(retried['executed'])
    prep_time=sum(d.segments[('P',h,0)]['duration'] for h in range(2))
    assert r4['duration_seconds']-r0['duration_seconds']==pytest.approx(3*prep_time)
    assert r0['location_counts']['2q']==36 and r4['location_counts']['2q']==102
    assert r0['location_counts']['reset']==16 and r4['location_counts']['reset']==64
    for q in DATA:
        e0=r0['idle_exposure_per_data_qubit'][str(q)];e4=r4['idle_exposure_per_data_qubit'][str(q)]
        assert e4['seconds']-e0['seconds']==pytest.approx(3*prep_time)
        assert e0['seconds']==pytest.approx(r0['duration_seconds']-2*JONES_MURALI_2026.ms_duration)
        for h in range(2):
            ops=d.segments[('P',h,0)]['ops']
            assert sum(op[2] for op in ops if op[0]=='IDLE' and op[1]==q)==pytest.approx(d.segments[('P',h,0)]['duration'])
    assert r4['ion_count']==15


def test_rejected_attempt_does_not_reset_data_or_hide_its_fault():
    # The Z detector has already run; a new Z on waiting data remains visible
    # at the endpoint while the following X-error detector leaves it alone.
    d=SteaneExtraction(IdleIonNoise(JONES_MURALI_2026),half_order='Z_FIRST')
    key=('P',1,0);ops=d.segments[key]['ops']
    idle=next(i for i,o in enumerate(ops) if o[0]=='IDLE' and o[1]==0)
    verify=next(i for i,o in enumerate(ops) if o[0] in ('M','MX'))
    faults={key:{idle:'Z',verify:'FLIP'}}
    out,f=d.run(faults=faults);_,t=d.run(faults=faults,tableau=True)
    assert f==t and out==(0,1) and not f['aborted']
    assert ('P',1,1) in f['executed'] and len(f['verification'])==3


@pytest.mark.parametrize('order',HALF_ORDERS)
def test_physical_four_rejections_abort_and_stop_costs(order):
    d=SteaneExtraction(scalar.UniformNoise(.001),half_order=order)
    faults={}
    for a in range(4):
        key=('P',0,a);idx=next(i for i,o in enumerate(d.segments[key]['ops']) if o[0] in ('M','MX'))
        faults[key]={idx:'FLIP'}
    for basis in ('0','+'):
        _,f=d.run(faults=faults,basis=basis);_,t=d.run(faults=faults,basis=basis,tableau=True)
        assert f==t and f['aborted'] and f['logical_failure']
        assert f['executed']==[('P',0,a) for a in range(4)]
        cost=d.resources(f['executed'])
        assert cost['location_counts']=={'reset':32,'2q':44,'measurement':4}


def test_random_codeword_readout_compares_syndrome_not_bits():
    # Published code span: XORing a codeword leaves H*m unchanged.
    clean=[0]*7
    for codeword in (0b0110011,0b1010101,0b1111000):
        bits=[(codeword>>q)&1 for q in DATA]
        assert syndrome(bits)==syndrome(clean)==0
        for q in DATA:
            noisy=bits.copy();noisy[q]^=1
            reference=clean.copy();reference[q]=1
            assert syndrome(noisy)==syndrome(reference)


def test_timing_free_H2_empty_transport_and_frozen_approved_results():
    d=SteaneExtraction(IonNoise(H2_LIKE));_,log=d.run()
    assert d.resources(log['executed'])['duration_seconds'] is None
    assert not any(op[0]=='IDLE' for _,_,op in d.locations)
    assert frozen_checks()=={'artifact_count':225,'all_byte_identical':True,'mismatches':[]}


@pytest.mark.parametrize('order',HALF_ORDERS)
@pytest.mark.parametrize('policy',['SEQ','PAR_MR'])
def test_memory_adapter_matches_scalar_with_retries_and_absorbing_aborts(order,policy):
    from steane_memory import SteaneMemoryRound,endpoint_score
    from memory_mc import memory_trial
    d=SteaneExtraction(IdleIonNoise(JONES_MURALI_2026),half_order=order,timing_policy=policy)
    a=SteaneMemoryRound(d);faults=[]
    for loc in d.locations:
        key,i,op=loc
        for payload,p in d.channel(loc).items():
            if p>0:faults.append({key:{i:payload}})
    # An actually reached retry fault after an explicit earlier rejection.
    first=('P',0,0);i=next(i for i,o in enumerate(d.segments[first]['ops']) if o[0] in ('M','MX'))
    retry=('P',0,1);ri=next(i for i,o in enumerate(d.segments[retry]['ops']) if o[0]=='IDLE' and o[1]==0)
    faults.append({first:{i:'FLIP'},retry:{ri:'Z'}})
    abort={}
    for k in range(4):
        key=('P',0,k);idx=next(i for i,o in enumerate(d.segments[key]['ops']) if o[0] in ('M','MX'))
        abort[key]={idx:'FLIP'}
    faults.append(abort)
    state=a.initial_state(len(faults));step=a.step(state,np.random.default_rng(8),faults=a.explicit_faults(faults))
    for j,f in enumerate(faults):
        out,log=d.run(faults=f)
        assert out==(int(step.state[0][j]),int(step.state[1][j]))
        assert log['aborted']==bool(step.aborted[j])
    assert step.aborted[-1] and not step.state[2][-1]
    follow=a.step(step.state,np.random.default_rng(8),faults={})
    assert not follow.state[2][-1] and follow.diagnostics['segment_execution_counts'][str(('P',0,0))]==len(faults)-1
    trial=memory_trial(a,1,4,np.random.default_rng(10),score=endpoint_score,
        fault_rounds=[a.explicit_faults([abort]),{}, {}, {}])
    assert trial['aborted']==trial['fails_X']==trial['fails_Z']==trial['fails']==1
    assert trial['round_totals']['attempts']==4
    assert all(v['attempts']==0 for v in trial['per_round'][1:])
