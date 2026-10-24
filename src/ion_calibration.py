"""Log-log MC calibration and joint binomial uncertainty, without extrapolation."""

import math
import numpy as np


def monotone_log_rates(rates, shots):
    """Weighted PAVA in log rate; weights approximate inverse log-rate variance.

    Point estimates usually need no adjustment. Bootstrap draws can contain
    small adjacent reversals; pooling enforces a single-valued inverse curve.
    """
    rates = np.asarray(rates, dtype=float)
    shots = np.asarray(shots, dtype=float)
    if np.any(~np.isfinite(rates)) or np.any(rates <= 0) or np.any(rates >= 1):
        raise ValueError('Calibration needs logical rates strictly between 0 and 1')
    weights = shots * rates / (1 - rates)
    if np.all(np.diff(rates) >= 0):
        return rates.copy()
    values = np.log(rates)
    blocks = []
    for i, (value, weight) in enumerate(zip(values, weights)):
        blocks.append([i, i + 1, value * weight, weight])
        while len(blocks) > 1 and blocks[-2][2] / blocks[-2][3] > blocks[-1][2] / blocks[-1][3]:
            right = blocks.pop()
            left = blocks.pop()
            blocks.append([left[0], right[1], left[2] + right[2], left[3] + right[3]])
    result = np.empty(len(values))
    for start, stop, total, weight in blocks:
        result[start:stop] = total / weight
    return np.exp(result)


def log_log_inverse(p, rates, target):
    p, rates = np.asarray(p, float), np.asarray(rates, float)
    target = np.asarray(target, float)
    if (len(p) < 2 or len(p) != len(rates) or np.any(np.diff(p) <= 0)
            or np.any(np.diff(rates) < 0) or np.any(p <= 0) or np.any(rates <= 0)
            or np.any(~np.isfinite(p)) or np.any(~np.isfinite(rates))):
        raise ValueError('Expected increasing positive p and nondecreasing positive rates')
    if np.any(~np.isfinite(target)) or np.any(target < rates[0]) or np.any(target > rates[-1]):
        raise ValueError('Target is not bracketed; extend the calibration grid')
    return np.exp(np.interp(np.log(target), np.log(rates), np.log(p)))


def pseudothreshold(p, rates):
    """First upward crossing of log(p_L)-log(p), interpolated in log p."""
    p, rates = np.asarray(p, float), np.asarray(rates, float)
    delta = np.log(rates) - np.log(p)
    for i in range(len(p) - 1):
        if delta[i] <= 0 <= delta[i + 1] and delta[i + 1] > delta[i]:
            fraction = -delta[i] / (delta[i + 1] - delta[i])
            return float(math.exp(math.log(p[i]) + fraction * math.log(p[i + 1] / p[i])))
    raise ValueError('Uniform MC pseudothreshold is not bracketed')


def _interval(samples):
    lo, hi = np.quantile(samples, [0.025, 0.975])
    return [float(lo), float(hi)]


def calibrate_targets(grid, targets, *, bootstrap_draws=20000, seed=302026):
    """Independent binomial resamples of target and grid, shared across ratios.

    targets maps a label to {'monte_carlo': {...}, 'p_2q': rate}. The same
    resampled grid yields p_eff and p_threshold in each draw, preserving their
    calibration covariance in p_eff/p_threshold confidence intervals.
    """
    if len(grid) < 6 or bootstrap_draws < 100:
        raise ValueError('Use at least six grid points and 100 bootstrap draws')
    grid = sorted(grid, key=lambda row: row['p'])
    p = np.array([row['p'] for row in grid])
    n = np.array([row['monte_carlo']['shots'] for row in grid])
    raw = np.array([row['monte_carlo']['fails'] / row['monte_carlo']['shots'] for row in grid])
    fitted = monotone_log_rates(raw, n)
    threshold = pseudothreshold(p, fitted)
    for target in targets.values():
        log_log_inverse(p, fitted, target['monte_carlo']['p_L'])
    rng = np.random.default_rng(seed)
    sampled_rates = rng.binomial(n, raw, size=(bootstrap_draws, len(p))) / n
    threshold_draws = np.full(bootstrap_draws, np.nan)
    adjusted = 0
    for i, row in enumerate(sampled_rates):
        try:
            fitted_row = monotone_log_rates(row, n)
            adjusted += int(np.any(np.diff(row) < 0))
            sampled_rates[i] = fitted_row
            threshold_draws[i] = pseudothreshold(p, fitted_row)
        except ValueError:
            continue
    valid_threshold = np.isfinite(threshold_draws)
    if np.mean(~valid_threshold) > 0.005:
        raise ValueError('Too many resampled thresholds unbracketed; extend the grid')
    output = {
        'method': 'Piecewise linear interpolation of log(p_L) versus log(p).',
        'confidence_method': '95% percentile parametric binomial bootstrap of target and uniform MC; shared grid draws for ratios.',
        'bootstrap_draws': bootstrap_draws, 'bootstrap_seed': seed,
        'monotone_adjustment': 'Inverse-log-variance-weighted PAVA on log rates if needed; flats use the rightmost endpoint for exact plateau values.',
        'point_grid_adjusted': bool(np.any(np.diff(raw) < 0)),
        'bootstrap_monotone_adjustment_fraction': adjusted / bootstrap_draws,
        'pseudothreshold': {'p': threshold, 'ci95': _interval(threshold_draws[valid_threshold]),
                           'valid_bootstrap_draws': int(sum(valid_threshold))},
        'grid_logical_rates_used': fitted.tolist(),
        'targets': {},
    }
    for label, target in targets.items():
        mc = target['monte_carlo']
        peff = float(log_log_inverse(p, fitted, mc['p_L']))
        target_draws = rng.binomial(mc['shots'], mc['p_L'], size=bootstrap_draws) / mc['shots']
        draws = np.full(bootstrap_draws, np.nan)
        for i, rate in enumerate(target_draws):
            try:
                draws[i] = log_log_inverse(p, sampled_rates[i], rate)
            except ValueError:
                continue
        valid = np.isfinite(draws) & valid_threshold
        if np.mean(~valid) > 0.005:
            raise ValueError(f'Target {label}: confidence interval not bracketed; extend grid')
        ci = _interval(draws[valid])
        neighborhood_start = min(max(int(np.searchsorted(p, peff)) - 3, 0), len(p) - 6)
        neighborhood = list(range(neighborhood_start, neighborhood_start + 6))
        output['targets'][label] = {
            'p_eff_MC': peff, 'ci95': ci,
            'p_eff_over_p_2q': peff / target['p_2q'],
            'p_eff_over_p_2q_ci95': [x / target['p_2q'] for x in ci],
            'p_eff_over_pseudothreshold': peff / threshold,
            'p_eff_over_pseudothreshold_ci95': _interval(draws[valid] / threshold_draws[valid]),
            'valid_bootstrap_draws': int(sum(valid)),
            'bracketing_six_point_grid': [float(p[i]) for i in neighborhood],
            'bracketing_six_point_shots': [int(n[i]) for i in neighborhood],
            'target_MC_shots': mc['shots'],
            'target_MC_fails': mc['fails'],
        }
    return output
