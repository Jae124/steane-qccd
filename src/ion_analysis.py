"""Stage 3: MC calibration, readout sensitivity, confidence intervals and figures.

The two expensive parts can be run independently with --part calibration and
--part measurements. --part report assembles validated cached data; default all
resumes matching caches. Every uniform point uses at least one million shots.
"""

import argparse
from dataclasses import asdict, replace
import hashlib
import inspect
import json
import math
from pathlib import Path
import random

from ion_calibration import calibrate_targets, log_log_inverse
from ion_noise import IdleIonNoise, IonNoise, check_channel_faults, run_preset
from ion_params import (JONES_MURALI_2026, MEASUREMENT_DURATION_SWEEP, PRESETS,
                        UNIFORM_CALIBRATION_P)
from ion_timing import idle_exposure
from steane_noise import UniformNoise, monte_carlo

ROOT = Path(__file__).resolve().parents[1]


def fingerprint():
    sources = [(ROOT / 'src' / name).read_bytes() for name in
               ('steane_ft.py', 'steane_noise.py', 'ion_timing.py', 'ion_params.py')]
    sources += [inspect.getsource(cls).encode() for cls in (IonNoise, IdleIonNoise)]
    return hashlib.sha256(b'\n'.join(sources)).hexdigest()


def save_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def _cache(path, configuration):
    if not path.exists():
        return {'configuration': configuration, 'points': {}}
    data = json.loads(path.read_text())
    if data['configuration'] != configuration:
        raise ValueError(f'{path}: cached simulation/configuration differs; choose a new cache directory')
    return data


def uniform_grid(cache_dir, shots, seed):
    if shots < 1_000_000:
        raise ValueError('Stage 3 uniform calibration requires >=1M shots at every point')
    config = {'simulation_fingerprint': fingerprint(), 'shots_per_point': shots,
              'seed_base': seed, 'p_grid': list(UNIFORM_CALIBRATION_P)}
    path = cache_dir / 'uniform_mc_calibration.json'
    data = _cache(path, config)
    if 'fault_checks' not in data:
        data['fault_checks'] = check_channel_faults(UniformNoise(UNIFORM_CALIBRATION_P[0]))
    checks = data['fault_checks']
    if checks['single_fault_failures'] or checks['mismatches']:
        save_json(path, data)
        raise RuntimeError('Uniform fault-tolerance/cross-check failure; stop for review')
    # Approved exact coefficients; uniform exhaustive regression remains in pytest.
    baseline_path = ROOT / 'results/stage0_channels.json'
    counts = json.loads(baseline_path.read_text())
    data['pair_coefficient_source_sha256'] = hashlib.sha256(baseline_path.read_bytes()).hexdigest()
    for i, p in enumerate(UNIFORM_CALIBRATION_P):
        key = str(p)
        if key in data['points']:
            continue
        point_seed = seed + 10000 + i
        print(f'Uniform grid: p={p:.8g}, {shots} shots, seed={point_seed}', flush=True)
        mc = monte_carlo(UniformNoise(p), shots, random.Random(point_seed), return_details=True)
        mc['seed'] = point_seed
        pair = counts['A']['any'] * p * p
        corrected = pair * (1 - p)**(counts['n_locations'] - 2)
        data['points'][key] = {'p': p, 'monte_carlo': mc, 'P1': 0,
            'P2': pair, 'P2_survival_corrected': corrected,
            'P2_over_MC': pair / mc['p_L'] if mc['p_L'] else None,
            'P2_survival_over_MC': corrected / mc['p_L'] if mc['p_L'] else None}
        save_json(path, data)
        print(f'  p_L={mc["p_L"]:.8g} +/- {mc["stderr"]:.3g}', flush=True)
    for point in data['points'].values():
        scale = point['p']**2
        survival = (1-point['p'])**(counts['n_locations']-2)
        for suffix, coefficient in [('type_pair', 'A_by_type_pair'), ('stage_pair', 'A_by_stage_pair')]:
            point['P2_by_' + suffix] = {key: value*scale for key, value in counts[coefficient].items()}
            point['P2_survival_by_' + suffix] = {key: value*scale*survival for key, value in counts[coefficient].items()}
    save_json(path, data)
    return data


def measurement_sweep(cache_dir, shots, seed):
    config = {'simulation_fingerprint': fingerprint(), 'shots_per_point': shots,
              'seed_base': seed, 'nominal_parameters': asdict(JONES_MURALI_2026),
              'measurement_durations_seconds': list(MEASUREMENT_DURATION_SWEEP)}
    path = cache_dir / 'measurement_duration_sweep.json'
    data = _cache(path, config)
    for i, duration in enumerate(MEASUREMENT_DURATION_SWEEP):
        key = str(duration)
        if key in data['points']:
            if data['points'][key]['status'] != 'PASS':
                raise RuntimeError('Cached fault-tolerance/cross-check failure; stop for review')
            continue
        params = replace(JONES_MURALI_2026, measurement_duration=duration)
        print(f'Readout sensitivity: measurement duration={duration*1e6:g} us', flush=True)
        result = run_preset(params, shots, seed + 1000 + 100*i, idle=True)
        data['points'][key] = result
        save_json(path, data)
        if result['status'] != 'PASS':
            raise RuntimeError('Readout sweep fault-tolerance/cross-check failure; stop for review')
        result['idle_exposure'] = idle_exposure(params, result['monte_carlo'])
        save_json(path, data)
    return data


def split_pair_attribution(counts):
    """One credit per pair: mixed-type pairs split equally between their types."""
    totals = dict.fromkeys(('2q', 'reset', 'meas', 'idle'), 0.0)
    for pair, value in counts['P2_by_type_pair'].items():
        a, b = pair.split('*')
        totals[a] += value / 2
        totals[b] += value / 2
    return totals


def make_figures(report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    figure_dir = ROOT / 'figures'
    figure_dir.mkdir(exist_ok=True)
    presets = list(report['presets'])
    types = ('2q', 'reset', 'meas', 'idle')
    colors = ('#2868a6', '#dc8a29', '#5c9970', '#9467bd')
    labels = ('2q', 'Reset (includes RX rotation)', 'Measurement (includes MX rotation)', 'Idle Z')
    fig, ax = plt.subplots(figsize=(9.0, 4.9), layout='constrained')
    positions = np.arange(len(presets))
    for i, (kind, color, label) in enumerate(zip(types, colors, labels)):
        values = [split_pair_attribution(report['presets'][p]['counts'])[kind] for p in presets]
        ax.bar(positions + (i-1.5)*0.18, values, width=0.17, color=color, label=label)
    ax.set_xticks(positions, [p.replace('_', ' ') + ('\nStage 1 only' if p == 'H2_LIKE' else '\nwith nominal idle Z')
                              for p in presets])
    ax.set_ylabel('Contribution to raw P2 (secondary attribution)')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.set_title('Steane EC pair attribution\nMixed-type pairs split equally; Monte Carlo is primary')
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.2)
    fig.savefig(figure_dir / 'ion_noise_breakdown.png', dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 5.2), layout='constrained')
    grid = report['uniform_calibration']['points']
    x = np.array([point['p'] for point in grid])
    y = np.array([point['monte_carlo']['p_L'] for point in grid])
    err = np.array([point['monte_carlo']['stderr'] for point in grid])
    ax.errorbar(x, y, yerr=1.96*err, fmt='o-', color='#2868a6', capsize=3,
                label='Uniform MC grid (95% sampling bars)')
    ax.plot(x, x, '--', color='0.4', label='p_L = p')
    for name, color in zip(presets, ('#b34d44', '#4d9377')):
        result = report['presets'][name]
        pe = result['p_eff']
        mc = result['monte_carlo']
        ax.errorbar(pe['p_eff_MC'], mc['p_L'],
                    xerr=[[pe['p_eff_MC']-pe['ci95'][0]], [pe['ci95'][1]-pe['p_eff_MC']]],
                    yerr=1.96*mc['stderr'], fmt='s', color=color, capsize=4,
                    label=name.replace('_', ' ') + (' target (Stage 1 only)' if name == 'H2_LIKE' else ' target'))
    ax.set(xscale='log', yscale='log', xlabel='Uniform per-location p', ylabel='Logical failure probability p_L',
           title='Monte Carlo effective uniform rates\nHorizontal intervals include target and calibration sampling error')
    ax.grid(alpha=0.2, which='both')
    ax.legend(fontsize=8)
    fig.savefig(figure_dir / 'ion_noise_calibration.png', dpi=200)
    plt.close(fig)


def assemble_report(cache_dir, output, *, preset='ALL', bootstrap_draws=20000, seed=2026):
    source = ROOT / 'results/stage2_ion_noise.json'
    prior = json.loads(source.read_text())
    grid_data = json.loads((cache_dir / 'uniform_mc_calibration.json').read_text())
    if grid_data['configuration']['simulation_fingerprint'] != fingerprint():
        raise ValueError('Calibration cache simulation fingerprint differs')
    grid = sorted(grid_data['points'].values(), key=lambda row: row['p'])
    if len(grid) != len(UNIFORM_CALIBRATION_P) or min(p['monte_carlo']['shots'] for p in grid) < 1_000_000:
        raise ValueError('Incomplete uniform grid or fewer than 1M shots per point')
    report = {'stage': 3, 'primary_estimate': 'Monte Carlo',
              'simulation_fingerprint': fingerprint(),
              'stage2_reference_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'uniform_calibration': {'configuration': grid_data['configuration'], 'points': grid,
                                      'fault_checks': grid_data['fault_checks']},
              'presets': {}}
    targets = {}
    if preset in ('ALL', 'JONES_MURALI_2026'):
        sweep = json.loads((cache_dir / 'measurement_duration_sweep.json').read_text())
        if sweep['configuration']['simulation_fingerprint'] != fingerprint():
            raise ValueError('Readout cache simulation fingerprint differs')
        if set(sweep['points']) != {str(d) for d in MEASUREMENT_DURATION_SWEEP}:
            raise ValueError('Incomplete measurement-duration sweep')
        for duration in MEASUREMENT_DURATION_SWEEP:
            point = sweep['points'][str(duration)]
            if point['status'] != 'PASS':
                raise ValueError('Readout point fault check failed; stop for review')
            targets[f'readout_{duration*1e6:g}us'] = {
                'monte_carlo': point['monte_carlo'], 'p_2q': point['parameters']['p_2q']}
        nominal_key = str(JONES_MURALI_2026.measurement_duration)
        # JSON copy isolates subsequent report annotations from cached input data.
        nominal = json.loads(json.dumps(sweep['points'][nominal_key]))
        nominal.update(stage=3, timing_applied=True,
                       T2_sweep=prior['presets']['JONES_MURALI_2026']['T2_sweep'],
                       T2_sweep_provenance='Approved Stage 2 reference, not newly calibrated T2 targets.')
        report['presets']['JONES_MURALI_2026'] = nominal
        report['measurement_duration_sweep'] = sweep
    if preset in ('ALL', 'H2_LIKE'):
        h2 = json.loads(json.dumps(prior['presets']['H2_LIKE']))
        h2.update(analysis_stage=3, scope='Stage 1 only', timing_applied=False,
                  primary_MC_provenance='Reused approved Stage 2 report timing-free H2 MC, 3M shots.')
        report['presets']['H2_LIKE'] = h2
        targets['H2_LIKE'] = {'monte_carlo': h2['monte_carlo'], 'p_2q': h2['parameters']['p_2q']}
    calibration = calibrate_targets(grid, targets, bootstrap_draws=bootstrap_draws, seed=seed+300000)
    p_values = [point['p'] for point in grid]
    report['uniform_round_trip_max_absolute_error'] = max(abs(
        float(log_log_inverse(p_values, calibration['grid_logical_rates_used'], rate))-p)
        for p, rate in zip(p_values, calibration['grid_logical_rates_used']))
    report['calibration_statistics'] = calibration
    report['uniform_pseudothreshold'] = calibration['pseudothreshold']
    A = json.loads((ROOT / 'results/stage0_channels.json').read_text())['A']['any']
    for name, result in report['presets'].items():
        key = (f'readout_{JONES_MURALI_2026.measurement_duration*1e6:g}us'
               if name == 'JONES_MURALI_2026' else 'H2_LIKE')
        result['p_eff'] = calibration['targets'][key]
        secondary = math.sqrt(result['counts']['P2'] / A)
        ci = result['p_eff']['ci95']
        result['p_eff_P2_secondary'] = {'p_eff': secondary, 'A_uniform': A,
            'p_eff_over_p_2q': secondary / result['parameters']['p_2q'],
            'within_MC_ci95': ci[0] <= secondary <= ci[1],
            'relative_to_MC': secondary / result['p_eff']['p_eff_MC']}
        result['P2_attribution_by_location_type'] = split_pair_attribution(result['counts'])
    if 'measurement_duration_sweep' in report:
        for duration, point in report['measurement_duration_sweep']['points'].items():
            point['p_eff'] = calibration['targets'][f'readout_{float(duration)*1e6:g}us']
            point['p_eff_P2_secondary'] = math.sqrt(point['counts']['P2'] / A)
        base = report['measurement_duration_sweep']['points']
        zero, nominal = base['0.0']['monte_carlo'], base[str(JONES_MURALI_2026.measurement_duration)]['monte_carlo']
        gate_only = prior['presets']['JONES_MURALI_2026']['stage1_reference']['monte_carlo']
        import numpy as np
        rng = np.random.default_rng(seed+400000)
        nominal_draw = rng.binomial(nominal['shots'], nominal['p_L'], bootstrap_draws)/nominal['shots']
        zero_draw = rng.binomial(zero['shots'], zero['p_L'], bootstrap_draws)/zero['shots']
        gate_draw = rng.binomial(gate_only['shots'], gate_only['p_L'], bootstrap_draws)/gate_only['shots']
        delta = nominal_draw-zero_draw
        valid = nominal_draw > gate_draw
        report['readout_sensitivity_comparison'] = {
            'p_L_400us_minus_0us': nominal['p_L'] - zero['p_L'],
            'difference_stderr': math.hypot(nominal['stderr'], zero['stderr']),
            'difference_ci95': np.quantile(delta, [0.025,0.975]).tolist(),
            'ratio_of_readout_increment_to_total_idle_increment': (nominal['p_L']-zero['p_L'])/(nominal['p_L']-gate_only['p_L']),
            'ratio_of_increments_ci95': (np.quantile(delta[valid]/(nominal_draw[valid]-gate_draw[valid]), [0.025,0.975]).tolist()
                                        if np.mean(valid) > 0.995 else None),
            'fraction_of_400us_p_L_removed_at_0us': 1 - zero['p_L'] / nominal['p_L'],
            'Stage1_no_idle_reference': gate_only,
            'note': 'Duration alone changes; measurement error remains nominal. Waiting of all live ions changes, including ancillas; this is not a data-only intervention.',
        }
    save_json(output, report)
    make_figures(report)
    return report


def run_stage3(*, preset='ALL', shots=3_000_000, grid_shots=1_000_000,
               bootstrap_draws=20000, seed=2026, output=None, cache_dir=None, part='all'):
    if shots <= 0 or bootstrap_draws < 100:
        raise ValueError('Positive shots and at least 100 bootstrap draws are required')
    cache_dir = ROOT / 'results' if cache_dir is None else Path(cache_dir)
    output = ROOT / 'results/ion_noise_results.json' if output is None else Path(output)
    if part in ('all', 'calibration'):
        uniform_grid(cache_dir, grid_shots, seed)
    if part in ('all', 'measurements') and (preset != 'H2_LIKE' or part == 'measurements'):
        measurement_sweep(cache_dir, shots, seed)
    if part in ('all', 'report'):
        report = assemble_report(cache_dir, output, preset=preset, bootstrap_draws=bootstrap_draws, seed=seed)
        print(json.dumps({'pseudothreshold': report['uniform_pseudothreshold'],
            'p_eff': {name: value['p_eff'] for name, value in report['presets'].items()}}, indent=2), flush=True)
        return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--part', choices=('all', 'calibration', 'measurements', 'report'), default='all')
    parser.add_argument('--shots', type=int, default=3_000_000)
    parser.add_argument('--grid-shots', type=int, default=1_000_000)
    parser.add_argument('--bootstrap-draws', type=int, default=20000)
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cache-dir', type=Path)
    args = parser.parse_args()
    run_stage3(**vars(args))
