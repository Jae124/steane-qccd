from dataclasses import replace
import itertools
import math

import pytest
import stim

from ion_noise import IonNoise, _state_signature, check_channel_faults
from ion_params import PRESETS, H2_LIKE, fidelity_to_depolarizing
from steane_ft import (DATA, apply_pauli, fault_list, prepare_logical, run_ops)
from steane_noise import (CIRCUITS, all_locations, count_faults,
                          ec_round_frame, logical_failure, validated_channel)


@pytest.mark.parametrize('n,fidelity,expected', [(1, 1, 0), (2, 1, 0),
                                               (1, 0.99, 0.015), (2, 0.99, 0.0125)])
def test_fidelity_conversion(n, fidelity, expected):
    assert fidelity_to_depolarizing(fidelity, n) == pytest.approx(expected)


@pytest.mark.parametrize('fidelity,n', [(1.1, 1), (-0.1, 2), (math.nan, 1),
                                       (0, 1), (0.99, 3)])
def test_fidelity_conversion_rejects_invalid(fidelity, n):
    with pytest.raises(ValueError):
        fidelity_to_depolarizing(fidelity, n)


@pytest.mark.parametrize('params', PRESETS.values(), ids=PRESETS.keys())
def test_ion_channels_and_faults(params):
    noise = IonNoise(params)
    for loc in all_locations():
        ch = validated_channel(noise, loc)
        assert 0 <= sum(ch.values()) <= 1
        assert all(prob >= 0 for prob in ch.values())
        ops = CIRCUITS[loc[:2]]
        assert all((loc[2], payload) in fault_list(ops) for payload in ch)
    checks = check_channel_faults(noise)
    assert checks['single_fault_runs'] == 1908
    assert checks['flagged_fault_runs'] == 1152
    assert checks['incoming_runs'] == 42
    assert checks['bare_triggered_pair_cross_check_runs'] == 756
    assert checks['single_fault_failures'] == 0, checks['failure_details']
    assert checks['mismatches'] == 0, checks['mismatch_details']


def test_rx_exact_pauli_composition():
    params = replace(H2_LIKE, p_reset=0.27, p_1q=0.39)
    # Independently enumerate reset Z followed by I/X/Y/Z rotation noise.
    paulis = {'I': (0, 0), 'X': (1, 0), 'Y': (1, 1), 'Z': (0, 1)}
    inverse = {v: k for k, v in paulis.items()}
    expected = dict.fromkeys(paulis, 0.0)
    for reset, pr in [('I', 1 - params.p_reset), ('Z', params.p_reset)]:
        for rotation, pg in [('I', 1 - params.p_1q)] + [(p, params.p_1q / 3) for p in 'XYZ']:
            r, g = paulis[reset], paulis[rotation]
            expected[inverse[(r[0] ^ g[0], r[1] ^ g[1])]] += pr * pg
    ch = IonNoise(params).channel(('F', 0, 0, ('RX', 8)))
    assert ch == pytest.approx({p: expected[p] for p in 'XYZ'})
    assert 1 - sum(ch.values()) == pytest.approx(expected['I'])


def test_mx_effective_flip_probability():
    params = replace(H2_LIKE, p_meas=0.2, p_1q=0.3)
    expected = 0.0
    for p, prob in [('I', 1 - params.p_1q)] + [(p, params.p_1q / 3) for p in 'XYZ']:
        for flipped, weight in [(False, 1 - params.p_meas), (True, params.p_meas)]:
            if (p in 'YZ') != flipped:
                expected += prob * weight
    assert IonNoise(params).channel(('F', 0, 0, ('MX', 8, 'flag'))) == pytest.approx({'FLIP': expected})


def test_terminal_mx_reduction_on_entangled_check_circuits():
    # Direct premeasurement Pauli vs a reported-bit flip, on both logical bases
    # with/without an incoming error. Compare both outputs and reduced data state.
    for ops in CIRCUITS.values():
        idx = next((i for i, op in enumerate(ops) if op[0] == 'MX'), None)
        if idx is None:  # Bare Z checks contain only a Z-basis measurement.
            continue
        for basis, incoming, pauli, meas_flip in itertools.product(('0', '+'), (False, True), 'IXYZ', (False, True)):
            outputs = []
            for explicit in (True, False):
                sim = stim.TableauSimulator()
                prepare_logical(sim, basis)
                if incoming:
                    apply_pauli(sim, [DATA[0]], 'Y')
                results = run_ops(sim, ops[:idx])
                if explicit and pauli != 'I':
                    apply_pauli(sim, [ops[idx][1]], pauli)
                flip = meas_flip if explicit else (meas_flip != (pauli in 'YZ'))
                results.update(run_ops(sim, ops[idx:], fault=(0, 'FLIP') if flip else None))
                outputs.append((results, _state_signature(sim, basis)))
            assert outputs[0] == outputs[1]


@pytest.mark.parametrize('spectator_rate', [0.0, 0.17, 1.0])
def test_survival_pair_product_with_nonuniform_spectator(spectator_rate):
    locs = all_locations()
    cx = [l for l in locs if l[0] == 'F' and l[3][0] == 'CX']
    selected = [cx[0], cx[2], locs[-1]]
    channels = [{'XI': 0.2, 'IZ': 0.05}, {'ZX': 0.3, 'YY': 0.1}, {'FLIP': spectator_rate}]
    class Sparse:
        def channel(self, loc):
            return channels[selected.index(loc)] if loc in selected else {}
    expected = raw = 0.0
    for a, b in itertools.combinations(range(3), 2):
        survival = math.prod(1 - sum(channels[i].values()) for i in range(3) if i not in (a, b))
        for va, pa in channels[a].items():
            for vb, pb in channels[b].items():
                faults = {}
                for loc, payload in [(selected[a], va), (selected[b], vb)]:
                    faults.setdefault(loc[:2], {})[loc[2]] = payload
                if any(logical_failure(*ec_round_frame(0, 0, faults))):
                    raw += pa * pb
                    expected += pa * pb * survival
    result = count_faults(Sparse())
    assert raw > 0
    assert 0 <= expected <= raw
    assert result['P2'] == pytest.approx(raw)
    assert result['P2_survival_corrected'] == pytest.approx(expected)
    assert sum(result['P2_survival_by_stage_pair'].values()) == pytest.approx(expected)
    assert sum(result['P2_survival_by_type_pair'].values()) == pytest.approx(expected)


@pytest.mark.slow
@pytest.mark.parametrize('params', PRESETS.values(), ids=PRESETS.keys())
def test_preset_pairs(params):
    result = count_faults(IonNoise(params))
    assert result['P1'] == 0
    assert result['single_fault_failures'] == 0
    assert 0 < result['P2_survival_corrected'] < result['P2']
    assert sum(result['P2_by_type_pair'].values()) == pytest.approx(result['P2'])
    assert sum(result['P2_survival_by_type_pair'].values()) == pytest.approx(result['P2_survival_corrected'])
