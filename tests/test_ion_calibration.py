from dataclasses import replace
import math

import numpy as np
import pytest

from ion_analysis import split_pair_attribution
from ion_calibration import (calibrate_targets, log_log_inverse,
                              monotone_log_rates, pseudothreshold)
from ion_params import JONES_MURALI_2026
from ion_timing import idle_exposure
from steane_ft import DATA


def synthetic_grid(shots=1_000_000):
    p = np.linspace(0.0008, 0.002, 8)
    rates = 750 * p**2
    return [{'p': float(x), 'monte_carlo': {'shots': shots, 'fails': round(y*shots),
                                         'p_L': round(y*shots)/shots}}
            for x, y in zip(p, rates)]


def test_exact_log_log_inverse_threshold_and_uniform_round_trip():
    p = np.geomspace(0.0007, 0.0023, 9)
    rates = 750*p*p
    for value in p:
        assert float(log_log_inverse(p, rates, 750*value*value)) == pytest.approx(value, rel=1e-14)
    assert float(log_log_inverse(p, rates, 750*0.00123**2)) == pytest.approx(0.00123, rel=1e-14)
    assert pseudothreshold(p, rates) == pytest.approx(1/750, rel=1e-14)
    with pytest.raises(ValueError, match='bracketed'):
        log_log_inverse(p, rates, rates[0]/2)
    with pytest.raises(ValueError, match='bracketed'):
        pseudothreshold(p, 100*p*p)


def test_monotone_pooling_only_when_needed():
    increasing = np.array([0.01,0.02,0.03])
    assert np.array_equal(monotone_log_rates(increasing, [1000]*3), increasing)
    fitted = monotone_log_rates([0.01,0.009,0.03], [1000]*3)
    assert np.all(np.diff(fitted) >= 0)
    assert fitted[0] == fitted[1]
    assert fitted[2] == pytest.approx(0.03)


def test_joint_binomial_uncertainty_includes_both_simulations():
    mc = {'shots': 3_000_000, 'fails': 5760, 'p_L': 0.00192}
    targets = {'target': {'monte_carlo': mc, 'p_2q': 0.001}}
    both = calibrate_targets(synthetic_grid(), targets, bootstrap_draws=2000, seed=8)
    only_grid = calibrate_targets(synthetic_grid(), {'target': {
        'monte_carlo': {'shots': 3_000_000_000, 'fails': 5_760_000, 'p_L': 0.00192},
        'p_2q': 0.001}}, bootstrap_draws=2000, seed=8)
    only_target = calibrate_targets(synthetic_grid(1_000_000_000), targets,
                                    bootstrap_draws=2000, seed=8)
    def width(report):
        lo,hi=report['targets']['target']['ci95']
        return hi-lo
    assert width(both) > width(only_grid)
    assert width(both) > width(only_target)
    item=both['targets']['target']
    assert item['ci95'][0] < item['p_eff_MC'] < item['ci95'][1]
    assert item['p_eff_MC'] == pytest.approx(0.0016, rel=0.001)
    assert item['p_eff_over_p_2q'] == pytest.approx(item['p_eff_MC']/0.001)
    assert item['p_eff_over_pseudothreshold'] == pytest.approx(
        item['p_eff_MC']/both['pseudothreshold']['p'])
    assert len(item['bracketing_six_point_grid']) == 6
    assert min(item['bracketing_six_point_shots']) >= 1_000_000


def test_bootstrap_reproducibility_and_no_extrapolation():
    targets={'x':{'monte_carlo':{'shots':1_000_000,'fails':1000,'p_L':0.001},'p_2q':0.001}}
    kwargs={'bootstrap_draws':200,'seed':19}
    a=calibrate_targets(synthetic_grid(),targets,**kwargs)
    assert a==calibrate_targets(synthetic_grid(),targets,**kwargs)
    with pytest.raises(ValueError):
        calibrate_targets(synthetic_grid()[:5],targets,**kwargs)
    targets['x']['monte_carlo']['p_L']=0.1
    with pytest.raises(ValueError,match='bracketed'):
        calibrate_targets(synthetic_grid(),targets,**kwargs)


def test_pair_type_attribution_conserves_probability():
    counts={'P2_by_type_pair':{'2q*2q':0.1,'2q*idle':0.2,'meas*reset':0.3}}
    split=split_pair_attribution(counts)
    assert split==pytest.approx({'2q':0.2,'idle':0.1,'meas':0.15,'reset':0.15})
    assert sum(split.values())==pytest.approx(sum(counts['P2_by_type_pair'].values()))


def test_adaptive_idle_exposure_uses_executed_path_histogram():
    hist={'shots':10,'triggers':4,'trigger_circuit_counts':[4,0,0,0,0,0], 'no_trigger_shots':6}
    result=idle_exposure(JONES_MURALI_2026,hist)
    assert sum(path['count'] for path in result['paths'].values())==10
    for q in DATA:
        metrics=result['monte_carlo_weighted_data'][str(q)]
        for name in ('seconds','summed_Z_probability'):
            expected=0.6*result['paths']['no_trigger']['data'][str(q)][name] + 0.4*result['paths']['trigger_F0']['data'][str(q)][name]
            assert metrics[name]==pytest.approx(expected)
            assert metrics[name+'_stderr']>0
    no_trigger={'shots':10,'triggers':0,'trigger_circuit_counts':[0]*6,'no_trigger_shots':10}
    result=idle_exposure(JONES_MURALI_2026,no_trigger)
    # All data are live for six flagged circuits; each data gate subtracts 40 us.
    assert result['paths']['no_trigger']['round_duration_seconds']==pytest.approx(0.0069)
    assert result['monte_carlo_weighted_data']['0']['seconds']==pytest.approx(0.00682)
    assert result['monte_carlo_weighted_data']['6']['seconds']==pytest.approx(0.00666)
    assert result['mean_readout_idle_seconds_per_data_qubit']==pytest.approx(0.0048)
    for q in DATA:
        assert result['monte_carlo_weighted_data'][str(q)]['seconds_stderr']==0
    infinite=idle_exposure(replace(JONES_MURALI_2026,T2=math.inf),no_trigger)
    assert infinite['monte_carlo_weighted_data']['0']['summed_Z_probability']==0
