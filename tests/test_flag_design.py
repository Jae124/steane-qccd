import hashlib
import json
import math
from pathlib import Path
import random

import pytest
import stim

import steane_ft as ft
import steane_noise as old
from check_order import frozen_artifacts, source_fingerprint
from flag_design import (CHECK_ORDERS, FlagDesign, TWO_QUBIT_ONLY, base_circuits,
                         ec_round, ec_round_frame, layer_timing, validate_round_order)
from ion_analysis import fingerprint
from ion_noise import IdleIonNoise, IonNoise, _state_signature
from ion_params import JONES_MURALI_2026, H2_LIKE
from ion_timing import operation_duration

ROOT=Path(__file__).parents[1]


def test_all_approved_artifacts_remain_byte_identical():
    assert frozen_artifacts()['all_byte_identical']
    assert frozen_artifacts()['artifact_count']==16
    assert fingerprint()=='481c4d8668af7a000bfed3b2e92163228375fd8a36d4af386dea49d5e8539dec'


@pytest.mark.parametrize('bad',[[],ft.ROUND[:-1],ft.ROUND+[('Z',0)], [('Z',0)]*6,
                               [('Z',True),*ft.ROUND[1:]],['Z','Z','Z','X','X','X']])
def test_invalid_round_permutations(bad):
    with pytest.raises(ValueError):validate_round_order(bad)


def test_default_order_delegates_without_altering_frozen_instrument():
    for model in [old.UniformNoise(0.001),IonNoise(H2_LIKE),IdleIonNoise(JONES_MURALI_2026)]:
        d=FlagDesign(model)
        explicit=FlagDesign(model,ft.ROUND)
        assert d.legacy and explicit.legacy
        assert d.circuits==getattr(model,'circuits',old.CIRCUITS)
        expected=old.monte_carlo(model,1000,random.Random(51),return_details=True)
        assert d.monte_carlo(1000,random.Random(51))==expected
        assert explicit.monte_carlo(1000,random.Random(51))==expected


def test_bare_syndrome_rows_are_independent_of_check_position():
    # Deliberately permute row order too, beyond A2's four named kind orders.
    order=tuple(reversed(ft.ROUND))
    for basis in ('0','+'):
        for q in ft.DATA:
            for p in 'XYZ':
                sim=stim.TableauSimulator();ft.prepare_logical(sim,basis)
                ft.apply_pauli(sim,[q],p)
                log=ec_round(sim,round_order=order)
                fx=(1<<q) if p in 'XY' else 0
                fz=(1<<q) if p in 'ZY' else 0
                trace={};after=ec_round_frame(fx,fz,{},round_order=order,log=trace)
                assert log==trace
                assert log['s_x']==(q+1 if p in 'XY' else 0)
                assert log['s_z']==(q+1 if p in 'ZY' else 0)
                assert old.logical_failure(*after)==(0,0)
                assert ft.output_is_correctable(sim,basis)


@pytest.mark.parametrize('name',CHECK_ORDERS)
def test_each_order_preserves_noiseless_state_and_program_gates(name):
    for policy in ('SEQ','PAR_MR'):
        d=FlagDesign(IdleIonNoise(JONES_MURALI_2026),CHECK_ORDERS[name],timing_policy=policy)
        for key,ops in d.circuits.items():
            assert [o for o in ops if o[0]!='IDLE']==base_circuits(CHECK_ORDERS[name])[key]
        for basis in ('0','+'):
            sim=stim.TableauSimulator();ft.prepare_logical(sim,basis)
            signature=_state_signature(sim,basis)
            log=ec_round(sim,circuits=d.circuits,round_order=d.round_order)
            assert log['trigger'] is None
            assert _state_signature(sim,basis)==signature
            assert ec_round_frame(0,0,{},circuits=d.circuits,round_order=d.round_order)==(0,0)


def test_parallel_liveness_time_conservation_and_bare_durations():
    p=JONES_MURALI_2026
    d=FlagDesign(IdleIonNoise(p),timing_policy='PAR_MR')
    for key,ops in d.base_circuits.items():
        starts,ends,total,_=layer_timing(ops,p,'PAR_MR')
        if key[0]=='F':
            assert total==pytest.approx(700e-6)
            assert starts[0]==starts[1]
            assert starts[len(ops)-2]==starts[len(ops)-1]
        else:
            assert total==pytest.approx(610e-6 if d.round_order[key[1]][0]=='Z' else 620e-6)
        for q in ft.DATA:
            active=sum(operation_duration(o,p) for o in ops if o[0]=='CX' and q in o[1:])
            wait=sum(o[2] for o in d.circuits[key] if o[0]=='IDLE' and o[1]==q)
            assert wait==pytest.approx(total-active)
    # The shorter R finishes before RX; its ion idles for the 5 us remainder.
    z=d.circuits[('F',0)]
    assert any(o[0]=='IDLE' and o[1]==7 and o[2]==pytest.approx(p.rotation_duration) for o in z)
    # Overlap has no effect on gate/reset/readout probabilities.
    for loc in old.all_locations(d.circuits):
        if loc[3][0]!='IDLE': assert d.channel(loc)==IonNoise(p).channel(loc)


def test_parallel_layers_do_not_merge_repeated_qubits_or_cross_transport():
    ops=[('R',7),('R',7),('RX',8),('M',7,'synd'),('MX',8,'flag')]
    starts,ends,_,_=layer_timing(ops,JONES_MURALI_2026,'PAR_MR')
    assert starts[1]>=ends[0]
    assert starts[1]==starts[2]
    starts,ends,_,_=layer_timing(ops,JONES_MURALI_2026,'PAR_MR',[(2,8,'shuttle')])
    assert starts[2]>=ends[1]+JONES_MURALI_2026.shuttle_duration-1e-15
    with pytest.raises(ValueError):FlagDesign(IonNoise(H2_LIKE),timing_policy='PAR_MR')


def test_transport_indices_and_exposure_follow_executed_order():
    order=CHECK_ORDERS['XZXZXZ']
    event=[('F',0,2,7,'shuttle')]
    for policy in ('SEQ','PAR_MR'):
        d=FlagDesign(IdleIonNoise(JONES_MURALI_2026),order,timing_policy=policy,
                     transport_events=event,heating_penalty={'shuttle':0.001})
        assert d.base_circuits[('F',0)][2]==('CX',7,3)
        idx=d.transport_plan.mapping[('F',0)][2]
        assert d.circuits[('F',0)][idx]==('CX',7,3)
        loc=next(l for l in old.all_locations(d.circuits) if l[:3]==('F',0,idx))
        assert sum(d.channel(loc).values())==pytest.approx(0.002)
        hist={'shots':1,'triggers':1,'trigger_circuit_counts':[1,0,0,0,0,0],'no_trigger_shots':0}
        exposure=d.exposure(hist)
        assert exposure['paths']['trigger_F0']['round_duration_seconds']==pytest.approx(
            sum(layer_timing(d.base_circuits[key],d.params,policy,d.transport_plan.grouped.get(key,()))[2]
                for key in [('F',0),*[('B',i) for i in range(6)]]))
    # A late data-ion event can be consumed in a later F check or in the
    # permuted bare round. Exercise the ordered conditional dispatch too.
    order=CHECK_ORDERS['XZXZXZ']
    late=[('F',1,len(base_circuits(order)[('F',1)]),0,'shuttle')]
    d=FlagDesign(IdleIonNoise(JONES_MURALI_2026),order,transport_events=late,
                 heating_penalty={'shuttle':0.001})
    assert d.dynamic_channels
    target=next(l for l in old.all_locations(d.circuits) if l[:2]==('B',4) and l[3][0]=='CX' and 0 in l[3][1:])
    assert sum(d.channel_for_trigger(target,1).values())==pytest.approx(0.002)
    assert sum(d.channel_for_trigger(target,4).values())==pytest.approx(0.001)
    trigger=next(l for l in old.all_locations(d.circuits) if l[:2]==('F',1) and l[3][0] in ('M','MX') and l[3][2]=='synd')
    for payload in old.physical_variants(target[3]):
        faults={trigger[:2]:{trigger[2]:'FLIP'},target[:2]:{target[2]:payload}}
        log={};after=ec_round_frame(0,0,faults,circuits=d.circuits,round_order=d.round_order,log=log)
        for basis in ('0','+'):
            sim=stim.TableauSimulator();ft.prepare_logical(sim,basis)
            trace=ec_round(sim,faults=faults,circuits=d.circuits,round_order=d.round_order)
            assert trace==log
            assert (not ft.output_is_correctable(sim,basis))==bool(old.logical_failure(*after)[basis=='+'])
    metrics=d.monte_carlo(1000,random.Random(9))
    assert sum(metrics['trigger_circuit_counts'])==metrics['triggers']
    # An explicit empty event override removes inherited transport completely.
    base=IdleIonNoise(JONES_MURALI_2026,event)
    cleared=FlagDesign(base,transport_events=[])
    assert cleared.circuits==IdleIonNoise(JONES_MURALI_2026).circuits


def test_two_qubit_only_has_exactly_the_requested_channel_support():
    model=TWO_QUBIT_ONLY(0.001)
    for loc in old.all_locations():
        assert sum(model.channel(loc).values())==pytest.approx(0.001 if loc[3][0]=='CX' else 0)


def test_exhaustive_matrix_covers_all_orders_presets_and_policies():
    data=json.loads((ROOT/'results/stageA1_fault_checks.json').read_text())
    assert data['source_fingerprint']==source_fingerprint()
    assert len(data['rows'])==36 and data['unique_exhaustive_checks']==20
    assert {r['order_name'] for r in data['rows']}==set(CHECK_ORDERS)
    for row in data['rows']:
        c=row['checks']
        assert c['single_fault_failures']==c['mismatches']==0
        assert c['incoming_runs']==42
        assert c['flagged_fault_runs']>0 and c['bare_only_fault_runs']>0


@pytest.mark.slow
def test_default_order_replays_frozen_stage3_and_stage4_empty_MC():
    d=FlagDesign(IdleIonNoise(JONES_MURALI_2026),round_order=ft.ROUND)
    actual=d.monte_carlo(3_000_000,random.Random(3326))
    actual['seed']=3326
    stage3=json.loads((ROOT/'results/ion_noise_results.json').read_text())['presets']['JONES_MURALI_2026']
    stage4=json.loads((ROOT/'results/stage4_transport_results.json').read_text())['presets']['JONES_MURALI_2026']
    assert actual==stage3['monte_carlo']==stage4['monte_carlo']
    assert actual['fails']==6121


@pytest.mark.slow
def test_default_order_pair_counts_are_frozen():
    expected=json.loads((ROOT/'results/ion_noise_results.json').read_text())['presets']['JONES_MURALI_2026']['counts']
    assert FlagDesign(IdleIonNoise(JONES_MURALI_2026),round_order=ft.ROUND).count_faults()==expected
