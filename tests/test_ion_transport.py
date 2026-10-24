from dataclasses import replace
import json
import math
from pathlib import Path
import random

import pytest

from ion_analysis import fingerprint
from ion_cache import sampling_manifest
from ion_noise import IdleIonNoise, check_channel_faults
from ion_params import JONES_MURALI_2026, H2_LIKE, SYNTHETIC_HEATING_PENALTIES
from ion_timing import idle_exposure, operation_duration
from ion_transport import (conditional_monte_carlo, synthetic_ancilla_events,
                           TransportPlan)
from steane_ft import DATA, A
from steane_noise import (CIRCUITS, all_locations, count_faults, ec_round_frame,
                          logical_failure, monte_carlo, physical_variants)

ROOT = Path(__file__).parents[1]


def test_empty_transport_exact_stage3_channels_and_sampling():
    base = IdleIonNoise(JONES_MURALI_2026)
    empty = IdleIonNoise(JONES_MURALI_2026, [], {'shuttle': 0.1})
    assert empty.circuits == base.circuits
    assert not empty.dynamic_channels
    for loc in all_locations(base.circuits):
        assert empty.channel(loc) == base.channel(loc)
    assert monte_carlo(empty, 10000, random.Random(3326), return_details=True) == monte_carlo(
        base, 10000, random.Random(3326), return_details=True)
    assert sampling_manifest() == '1edaef4078655b02fd76d68086b7d984487db915e8e7ecc1e889f6ccbb8b7bd5'
    grid = json.loads((ROOT/'results/uniform_mc_calibration.json').read_text())
    assert fingerprint() == grid['configuration']['simulation_fingerprint']


def test_all_event_kinds_charge_every_live_qubit_and_consumed_heating():
    params = JONES_MURALI_2026
    events = [('F', 0, 2, A, kind) for kind in ('shuttle','split','merge','junction')]
    penalties = {'shuttle':0.001,'split':0.002,'merge':0.003,'junction':0.004}
    n = IdleIonNoise(params, events, penalties)
    base = IdleIonNoise(params)
    delta = sum(getattr(params, kind+'_duration') for kind in penalties)
    for q in [*DATA, 7, 8]:
        def wait(model):
            return sum(op[2] for op in model.circuits[('F',0)] if op[0]=='IDLE' and op[1]==q)
        assert wait(n)-wait(base) == pytest.approx(delta)
    first = n.transport_plan.mapping[('F',0)][2]
    loc = next(l for l in all_locations(n.circuits) if l[:3]==('F',0,first))
    assert sum(n.channel(loc).values()) == pytest.approx(params.p_2q+sum(penalties.values()))
    other = [l for l in all_locations(n.circuits) if l[3][0]=='CX' and l!=loc]
    assert all(sum(n.channel(l).values()) == pytest.approx(params.p_2q) for l in other)


def test_live_windows_before_reset_after_measurement_and_boundary_carry():
    p=JONES_MURALI_2026
    # Before reset neither ancilla is live. After both measurements only data
    # remain live. Moving an unprepared ancilla still costs all live data time.
    events=[('F',0,0,A,'shuttle'),('F',0,len(CIRCUITS[('F',0)]),A,'junction')]
    n=IdleIonNoise(p,events,{'shuttle':0.001,'junction':0.002})
    base=IdleIonNoise(p)
    def wait(model,q): return sum(o[2] for o in model.circuits[('F',0)] if o[0]=='IDLE' and o[1]==q)
    for q in DATA:
        assert wait(n,q)-wait(base,q)==pytest.approx(p.shuttle_duration+p.junction_duration)
    for q in (7,8): assert wait(n,q)==wait(base,q)
    # End-of-F0 heating follows next executed CX: B0 on a trigger, F1 otherwise.
    for stage, circuit in [('F',1),('B',0)]:
        first = next(i for i,o in enumerate(CIRCUITS[(stage,circuit)]) if o[0]=='CX')
        idx=n.transport_plan.mapping[(stage,circuit)][first]
        loc=next(l for l in all_locations(n.circuits) if l[:3]==(stage,circuit,idx))
        assert sum(n.channel_for_trigger(loc,0).values())==pytest.approx(p.p_2q+0.002)


def test_heating_branch_selects_next_executed_gate_and_reset_does_not_cool():
    p=JONES_MURALI_2026
    event=('F',1,len(CIRCUITS[('F',1)]),4,'shuttle')
    n=IdleIonNoise(p,[event],{'shuttle':0.003})
    assert n.dynamic_channels
    original=next(i for i,o in enumerate(CIRCUITS[('B',0)]) if o[0]=='CX' and 4 in o[1:])
    idx=n.transport_plan.mapping[('B',0)][original]
    loc=next(l for l in all_locations(n.circuits) if l[:3]==('B',0,idx))
    assert sum(n.channel_for_trigger(loc,0).values())==pytest.approx(p.p_2q)
    assert sum(n.channel_for_trigger(loc,1).values())==pytest.approx(p.p_2q+0.003)
    assert sum(n.channel_for_trigger(loc,2).values())==pytest.approx(p.p_2q)
    # Reset of A in the next circuit does not consume its pending heating.
    a=IdleIonNoise(p,[('F',0,len(CIRCUITS[('F',0)]),A,'shuttle')],{'shuttle':0.003})
    first=next(l for l in all_locations(a.circuits) if l[:2]==('B',0) and l[3][0]=='CX')
    assert sum(a.channel_for_trigger(first,0).values())==pytest.approx(p.p_2q+0.003)


@pytest.mark.parametrize('event', [
    ('Q',0,0,7,'shuttle'),('F',6,0,7,'shuttle'),('F',0,-1,7,'shuttle'),
    ('F',0,100,7,'shuttle'),('F',0,0,9,'shuttle'),('F',0,0,7,'fly'),('F',0,0),
])
def test_invalid_events(event):
    with pytest.raises(ValueError): IdleIonNoise(JONES_MURALI_2026,[event])


@pytest.mark.parametrize('penalty', [{'shuttle':-0.1},{'split':math.nan},{'bad':0}, {'shuttle':1.0}])
def test_invalid_or_overflowing_heating(penalty):
    with pytest.raises(ValueError):
        IdleIonNoise(JONES_MURALI_2026,[('F',0,2,7,'shuttle')],penalty)
    with pytest.raises(ValueError): IdleIonNoise(H2_LIKE,[('F',0,1,7,'shuttle')])


def test_synthetic_single_faults_and_transport_exposure():
    n=IdleIonNoise(JONES_MURALI_2026,synthetic_ancilla_events(),{'shuttle':0.001})
    assert not n.dynamic_channels  # symmetric fixture has equal B contexts
    checks=check_channel_faults(n)
    assert checks['flagged_fault_runs']==1344
    assert checks['bare_only_fault_runs']==912
    assert checks['single_fault_failures']==checks['mismatches']==0, checks
    hist={'shots':10,'triggers':0,'trigger_circuit_counts':[0]*6,'no_trigger_shots':10}
    expo=idle_exposure(n.params,hist,circuits=n.circuits,extra_durations=n.transport_plan.extra_durations)
    base=idle_exposure(n.params,hist)
    for q in DATA:
        assert expo['monte_carlo_weighted_data'][str(q)]['seconds']-base['monte_carlo_weighted_data'][str(q)]['seconds']==pytest.approx(12*n.params.shuttle_duration)
    assert expo['paths']['no_trigger']['round_duration_seconds']-base['paths']['no_trigger']['round_duration_seconds']==pytest.approx(12*n.params.shuttle_duration)
    # Conditional sampler has identical RNG results when contexts are equal.
    assert conditional_monte_carlo(n,5000,random.Random(7),return_details=True)==monte_carlo(n,5000,random.Random(7),return_details=True)


def sparse_branch_model():
    p=replace(JONES_MURALI_2026,p_2q=0,p_1q=0,p_reset=0,p_meas=0,T2=math.inf)
    inner=IdleIonNoise(p,[('F',1,len(CIRCUITS[('F',1)]),4,'shuttle')],{'shuttle':1.0})
    trigger=next(l for l in all_locations(inner.circuits) if l[:2]==('F',1) and l[3][0] in ('M','MX') and l[3][2]=='synd')
    target=next(l for l in all_locations(inner.circuits) if l[:2]==('B',0) and l[3][0]=='CX' and 4 in l[3][1:])
    class Sparse:
        dynamic_channels=True
        circuits=inner.circuits
        def channel(self,loc): return self.channel_for_trigger(loc,None)
        def channel_for_trigger(self,loc,stop):
            if loc==trigger: return {'FLIP':1.0}
            return {v:x for v,x in inner.channel_for_trigger(loc,stop).items() if x>0}
    return Sparse(),trigger,target


def test_conditional_sampler_and_pairs_against_sparse_exact_oracle():
    n,trigger,target=sparse_branch_model()
    expected=0.0
    for v in physical_variants(target[3]):
        f={trigger[:2]:{trigger[2]:'FLIP'},target[:2]:{target[2]:v}}
        expected+=any(logical_failure(*ec_round_frame(0,0,f,circuits=n.circuits)))/15
    assert expected>0
    result=count_faults(n)
    assert result['P1']==0
    assert result['P2']==pytest.approx(expected)
    mc=monte_carlo(n,10000,random.Random(9),return_details=True)
    assert mc['trigger_circuit_counts']==[0,10000,0,0,0,0]
    assert mc['p_L']==pytest.approx(expected,abs=0.02)


@pytest.mark.slow
def test_empty_transport_counts_exact_stage3_and_synthetic_P2_monotone():
    old=json.loads((ROOT/'results/ion_noise_results.json').read_text())['presets']['JONES_MURALI_2026']['counts']
    assert count_faults(IdleIonNoise(JONES_MURALI_2026,[]))==old
    values=[]
    for penalty in SYNTHETIC_HEATING_PENALTIES:
        n=IdleIonNoise(JONES_MURALI_2026,synthetic_ancilla_events(),{'shuttle':penalty})
        result=count_faults(n)
        assert result['P1']==0
        values.append(result['P2'])
    assert all(a<b for a,b in zip(values,values[1:]))
