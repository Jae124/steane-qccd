from dataclasses import replace
import math
import random

import pytest
import stim

from ion_noise import (IdleIonNoise, IonNoise, _state_signature,
                       check_channel_faults)
from ion_params import H2_LIKE, JONES_MURALI_2026
from ion_timing import idle_probability, insert_idles, operation_duration
from steane_ft import DATA, ec_round, ops_to_stim, prepare_logical
from steane_noise import (CIRCUITS, all_locations, count_faults,
                          ec_round_frame, monte_carlo, trigger_contributors)


def test_timing_noiseless_and_total_waits():
    noise = IdleIonNoise(JONES_MURALI_2026)
    for key, ops in noise.circuits.items():
        base = CIRCUITS[key]
        assert [op for op in ops if op[0] != 'IDLE'] == base
        assert all(op[2] > 0 for op in ops if op[0] == 'IDLE')
        elapsed = math.fsum(operation_duration(op, noise.params) for op in base)
        for q in DATA:
            wait = sum(op[2] for op in ops if op[0] == 'IDLE' and op[1] == q)
            active = sum(operation_duration(op, noise.params) for op in base
                         if op[0] == 'CX' and q in op[1:])
            assert wait == pytest.approx(elapsed - active)
        assert isinstance(ops_to_stim(ops), stim.Circuit)
    for basis in ('0', '+'):
        sim = stim.TableauSimulator()
        prepare_logical(sim, basis)
        before = _state_signature(sim, basis)
        log = ec_round(sim, circuits=noise.circuits)
        assert log['trigger'] is None
        assert _state_signature(sim, basis) == before
        assert ec_round_frame(0, 0, {}, circuits=noise.circuits) == (0, 0)


def test_idle_probability_and_boundary_composition():
    p = JONES_MURALI_2026
    wait = p.measurement_duration
    prob = idle_probability(wait, p.T2)
    assert prob == pytest.approx(9.0900826947e-5, rel=1e-8)
    # Splitting an untouched cross-circuit gap at the decision boundary is exact.
    a, b = idle_probability(wait, p.T2), idle_probability(3 * wait, p.T2)
    assert a + b - 2 * a * b == pytest.approx(idle_probability(4 * wait, p.T2))
    assert idle_probability(wait, math.inf) == 0
    assert idle_probability(0, p.T2) == 0
    with pytest.raises(ValueError):
        idle_probability(-1, p.T2)
    with pytest.raises(ValueError):
        idle_probability(wait, None)
    with pytest.raises(ValueError):
        IdleIonNoise(H2_LIKE)


def test_ancilla_live_window_and_no_time_charged_for_skipped_circuits():
    params = JONES_MURALI_2026
    # A is not live before its reset; F's lifetime excludes the wait before
    # its reset and after its measurement. R/RX composite duration is included.
    ops = [('R', 7), ('RX', 8), ('CX', 8, 7), ('M', 7, 'synd'), ('MX', 8, 'flag')]
    timed = insert_idles(ops, params)
    a_wait = sum(op[2] for op in timed if op[0] == 'IDLE' and op[1] == 7)
    f_wait = sum(op[2] for op in timed if op[0] == 'IDLE' and op[1] == 8)
    assert a_wait == params.reset_duration + params.rotation_duration
    assert f_wait == params.measurement_duration
    noise = IdleIonNoise(params)
    # A single fault in any bare-only idle is never executed on noiseless input.
    for stage, circuit, idx, op in all_locations(noise.circuits):
        if stage == 'B' and op[0] == 'IDLE':
            log = {}
            assert ec_round_frame(0, 0, {(stage, circuit): {idx: 'Z'}},
                                  circuits=noise.circuits, log=log) == (0, 0)
            assert log['trigger'] is None


def test_all_idle_faults_and_cross_check():
    checks = check_channel_faults(IdleIonNoise(JONES_MURALI_2026))
    assert checks['flagged_fault_runs'] == 1332
    assert checks['bare_only_fault_runs'] == 888
    assert checks['flagged_idle_fault_runs'] == 180
    assert checks['bare_only_idle_fault_runs'] == 132
    assert checks['flagged_fault_failures'] == 0, checks['failure_details']
    assert checks['bare_only_fault_failures'] == 0, checks['failure_details']
    assert checks['incoming_failures'] == 0, checks['failure_details']
    assert checks['mismatches'] == 0, checks['mismatch_details']


def test_infinite_T2_exact_stage1_sampling_and_channels():
    params = replace(JONES_MURALI_2026, T2=math.inf)
    idle = IdleIonNoise(params)
    baseline = IonNoise(params)
    assert idle.circuits is CIRCUITS
    assert [idle.channel(l) for l in all_locations()] == [baseline.channel(l) for l in all_locations()]
    assert monte_carlo(idle, 10000, random.Random(2026), return_details=True) == monte_carlo(
        baseline, 10000, random.Random(2026), return_details=True)


def test_trigger_attribution_single_categories_and_mixed_cancellation():
    circuits = IdleIonNoise(JONES_MURALI_2026).circuits
    locs = all_locations(circuits)
    seen = set()
    for stage, circuit, idx, op in locs:
        if stage != 'F':
            continue
        # Search single faults rather than assuming a particular check index.
        from steane_ft import fault_list
        for _, payload in fault_list([op]):
            faults = {(stage, circuit): {idx: payload}}
            log = {}
            ec_round_frame(0, 0, faults, circuits=circuits, log=log)
            if log['trigger'] is not None:
                result = trigger_contributors(faults, circuits, log)
                assert len(result) == 1
                seen.update(result)
    assert seen == {'reset', 'measurement', '2q', 'idle'}
    sample = {('F', 0): {next(i for i, o in enumerate(circuits[('F', 0)])
                                 if o[0] == 'M'): 'FLIP'}}
    log = {}
    ec_round_frame(0, 0, sample, circuits=circuits, log=log)
    assert trigger_contributors(sample, circuits, log) == ['measurement']
    ops = circuits[('F', 0)]
    reset = next(i for i, o in enumerate(ops) if o[0] == 'R')
    cx = next(i for i, o in enumerate(ops) if o[0] == 'CX')
    # Reset X and CX target X cancel exactly, leaving no trigger.
    cancelled = {('F', 0): {reset: 'X', cx: 'IX'}}
    ec_round_frame(0, 0, cancelled, circuits=circuits, log=log)
    assert log['trigger'] is None
    # Adding readout FLIP makes XOR(1,1,1)=1: all three types contribute.
    triple = {('F', 0): {**cancelled[('F', 0)], **sample[('F', 0)]}}
    ec_round_frame(0, 0, triple, circuits=circuits, log=log)
    assert set(trigger_contributors(triple, circuits, log)) == {'reset', '2q', 'measurement'}
    # When the syndrome contributions cancel, only an independent flag idle
    # error contributes to the actual nonzero bit; don't credit the other types.
    flag_idle = max(i for i, o in enumerate(ops) if o[0] == 'IDLE' and o[1] == 8)
    idle_trigger = {('F', 0): {**cancelled[('F', 0)], flag_idle: 'Z'}}
    ec_round_frame(0, 0, idle_trigger, circuits=circuits, log=log)
    assert trigger_contributors(idle_trigger, circuits, log) == ['idle']

    class Deterministic:
        def __init__(self):
            self.circuits = circuits
        def channel(self, loc):
            payload = triple.get(loc[:2], {}).get(loc[2])
            return {payload: 1.0} if payload is not None else {}
    metrics = monte_carlo(Deterministic(), 10, random.Random(1), return_details=True)
    assert metrics['trigger_rate'] == 1
    assert metrics['trigger_cause_counts']['mixed'] == 10
    assert metrics['trigger_fractional_type_rates']['reset'] == pytest.approx(1/3)


def test_metric_accounting_preserves_seeded_failure_sampling():
    noise = IonNoise(H2_LIKE)
    result = monte_carlo(noise, 10000, random.Random(2026), return_details=True)
    assert result['fails'] == monte_carlo(noise, 10000, random.Random(2026))
    assert sum(result['trigger_cause_counts'].values()) == result['triggers']
    assert sum(result['trigger_fractional_type_rates'].values()) == pytest.approx(result['trigger_rate'])
    assert 0 <= result['failure_share_through_bare'] <= 1


@pytest.mark.slow
def test_infinite_T2_pair_count_exactly_matches_stage1():
    import json
    from pathlib import Path
    expected = json.loads((Path(__file__).parents[1] / 'results/stage1_ion_noise.json').read_text())
    expected = expected['presets']['JONES_MURALI_2026']['counts']
    actual = count_faults(IdleIonNoise(replace(JONES_MURALI_2026, T2=math.inf)))
    assert actual == expected


def test_lumped_Z_matches_injection_at_beginning_of_each_wait():
    from steane_noise import _apply_payload, run_ops_frame
    noise = IdleIonNoise(JONES_MURALI_2026)
    for key, ops in noise.circuits.items():
        base = CIRCUITS[key]
        for idx, op in enumerate(ops):
            if op[0] != 'IDLE':
                continue
            q = op[1]
            prefix = [o for o in ops[:idx] if o[0] != 'IDLE']
            last_touch = max((i for i, o in enumerate(prefix)
                              if q in (o[1:3] if o[0] == 'CX' else (o[1],))), default=-1)
            cut = last_touch + 1
            fx, fz, before = run_ops_frame(0, 0, base[:cut], {})
            fx, fz = _apply_payload(fx, fz, [q], 'Z')
            fx, fz, after = run_ops_frame(fx, fz, base[cut:], {})
            early = fx, fz, {**before, **after}
            late = run_ops_frame(0, 0, ops, {idx: 'Z'})
            assert early == late
