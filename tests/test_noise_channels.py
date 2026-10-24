import json
import random
from pathlib import Path

import pytest

from steane_ft import single_fault_test
from steane_noise import (NoiseModel, UniformNoise, all_locations, count_faults,
                          cross_check, monte_carlo, physical_variants,
                          validated_channel, ec_round_frame, logical_failure)

BASELINE = json.loads((Path(__file__).parents[1] / 'results/steane_ec_results.json').read_text())


def test_uniform_channels_and_legacy_api():
    for loc in all_locations():
        ch = validated_channel(UniformNoise(0.01), loc)
        assert all(p >= 0 for p in ch.values())
        assert sum(ch.values()) == pytest.approx(0.01)
        assert list(ch) == physical_variants(loc[3])
        assert ch == NoiseModel(0.01).channel(loc)
    assert len(all_locations()) == 96


@pytest.mark.parametrize('channel', [{'X': -0.1}, {'X': 1.1}, {'X': float('nan')}, {'BAD': 0.1}])
def test_invalid_channels(channel):
    class Invalid:
        def channel(self, loc):
            return channel
    with pytest.raises(ValueError):
        validated_channel(Invalid(), all_locations()[0])


def test_single_fault_and_cross_check():
    assert single_fault_test(True, verbose_failures=0) == 0
    assert single_fault_test(False, verbose_failures=0) == 48
    assert cross_check() == {'runs': 1176, 'mismatches': 0}


def test_nonuniform_channel_weighting():
    # A sparse model with unequal payload weights tests the generic interface
    # independently of the uniform-rate adapter.
    locs = all_locations()
    cnot_locs = [l for l in locs if l[0] == 'F' and l[3][0] == 'CX']
    selected = [cnot_locs[0], cnot_locs[2]]
    distributions = [{'XI': 0.2, 'IZ': 0.05}, {'ZX': 0.3, 'YY': 0.1}]
    class Sparse:
        def channel(self, loc):
            return distributions[selected.index(loc)] if loc in selected else {}
    expected = 0.0
    for va, pa in distributions[0].items():
        for vb, pb in distributions[1].items():
            faults = {}
            for loc, payload in zip(selected, [va, vb]):
                faults.setdefault(loc[:2], {})[loc[2]] = payload
            if any(logical_failure(*ec_round_frame(0, 0, faults))):
                expected += pa * pb
    assert expected > 0
    result = count_faults(Sparse())
    assert result['P1'] == 0
    assert result['P2'] == pytest.approx(expected)
    assert 'A' not in result
    # Exercise unequal-weight sampling and compare with its exact two-location result.
    shots = 30000
    measured = monte_carlo(Sparse(), shots, random.Random(12)) / shots
    assert measured == pytest.approx(expected, abs=0.01)


@pytest.mark.slow
def test_exhaustive_uniform_pairs():
    p = 5e-4
    result = count_faults(UniformNoise(p))
    baseline = BASELINE['counts']
    for key in ('n_locations', 'location_counts', 'single_fault_failures', 'malignant_pairs'):
        assert result[key] == baseline[key]
    for key in ('A', 'A_by_stage_pair', 'A_by_type_pair'):
        assert result[key] == baseline[key]
    assert result['P1'] == 0
    assert result['P2'] == pytest.approx(result['A']['any'] * p**2, rel=2e-12)
    assert sum(result['P2_by_type_pair'].values()) == pytest.approx(result['P2'])
    assert sum(result['P2_by_stage_pair'].values()) == pytest.approx(result['P2'])
    assert result['P2_survival_corrected'] == pytest.approx(
        result['P2'] * (1 - p)**94, rel=2e-12)
    assert result['P2_survival_mean_approx'] == pytest.approx(
        result['P2_survival_corrected'], rel=2e-12)


@pytest.mark.slow
def test_seeded_uniform_monte_carlo():
    assert monte_carlo(UniformNoise(5e-4), 3_000_000, random.Random(2026)) == 542
