"""Stage 1: operation-specific ion channels and checkpoint report.

Run from src/ or the repository root. No idle or transport noise is applied.
"""

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import random

import stim

from ion_params import PRESETS
from steane_ft import (DATA, SUPPORTS, apply_pauli, ec_round, fault_list,
                       output_is_correctable, pauli_on, prepare_logical)
from steane_noise import (DATA_MASK, SUPPORT_MASKS, _parity, all_locations,
                          count_faults, ec_round_frame, logical_failure, monte_carlo,
                          physical_variants, validated_channel)


class IonNoise:
    """Exact local Pauli composition; terminal MX rotation is marginalized."""

    def __init__(self, params):
        self.params = params

    def channel(self, loc):
        op = loc[3][0]
        p = self.params
        if op == 'CX':
            return {v: p.p_2q / 15 for v in physical_variants(loc[3])}
        if op == 'R':
            return {'X': p.p_reset}
        if op == 'RX':
            # Compose reset Z with independent rotation I/X/Y/Z, modulo phase.
            a = p.p_1q / 3
            return {'X': a, 'Y': a,
                    'Z': (1 - p.p_reset) * a + p.p_reset * (1 - p.p_1q)}
        if op == 'M':
            return {'FLIP': p.p_meas}
        if op == 'MX':
            # X commutes with MX; Y/Z flip it. Independent flips combine by XOR.
            rotation_flip = 2 * p.p_1q / 3
            return {'FLIP': p.p_meas + rotation_flip - 2 * p.p_meas * rotation_flip}
        raise ValueError(f'Unsupported operation {op}')


def _state_signature(sim, basis):
    """Stabilizer syndromes and logical parity BEFORE ideal output decoding."""
    return tuple(sim.peek_observable_expectation(pauli_on(s, pauli))
                 for pauli in ('Z', 'X') for s in SUPPORTS) + (
        sim.peek_observable_expectation(pauli_on(DATA, 'Z' if basis == '0' else 'X')),)


def _frame_signature(fx, fz, basis):
    return tuple(1 - 2 * _parity(frame & mask)
                 for frame in (fx, fz) for mask in SUPPORT_MASKS) + (
        1 - 2 * _parity((fx if basis == '0' else fz) & DATA_MASK),)


def check_channel_faults(noise):
    """Check every supported single fault in both stages, plus incoming errors.

    Bare-only faults are skipped by a noiseless flagged round. Separately force
    a flagged measurement trigger and check every bare fault against tableau;
    these two-fault trials are cross-checks, not single-fault FT claims.
    """
    failures, mismatches = [], []
    runs = incoming_runs = bare_pair_runs = 0
    flagged_runs = bare_runs = 0
    locs = all_locations()
    allowed = {}
    for loc in locs:
        allowed.setdefault(loc[:2], []).append(loc[3])
    allowed = {key: set(fault_list(ops)) for key, ops in allowed.items()}
    trigger_loc = next(l for l in locs
                       if l[:2] == ('F', 0) and l[3][0] == 'M')

    def trial(faults, loc, incoming=None, require_ft=True):
        nonlocal runs, incoming_runs, bare_pair_runs, flagged_runs, bare_runs
        fx = fz = 0
        if incoming:
            q, p = incoming
            fx = (1 << q) if p in 'XY' else 0
            fz = (1 << q) if p in 'YZ' else 0
        fx, fz = ec_round_frame(fx, fz, faults)
        frame_failure = logical_failure(fx, fz)
        for basis, frame_fail in (('0', frame_failure[0]), ('+', frame_failure[1])):
            sim = stim.TableauSimulator()
            prepare_logical(sim, basis)
            if incoming:
                apply_pauli(sim, [incoming[0]], incoming[1])
            log = ec_round(sim, faults=faults)
            signature = _state_signature(sim, basis)
            tableau_fail = not output_is_correctable(sim, basis)
            record = {'basis': basis, 'fault': loc, 'incoming': incoming,
                      'trigger': log['trigger'], 's_x': log.get('s_x'),
                      's_z': log.get('s_z'), 'x_fix': log.get('x_fix', []),
                      'z_fix': log.get('z_fix', []),
                      'residual_frame': {'X': fx & DATA_MASK, 'Z': fz & DATA_MASK},
                      'tableau_signature': signature,
                      'frame_signature': _frame_signature(fx, fz, basis)}
            if require_ft and (tableau_fail or frame_fail):
                failures.append(record)
            if signature != record['frame_signature'] or bool(frame_fail) != tableau_fail:
                mismatches.append(record)
            if incoming:
                incoming_runs += 1
            elif require_ft:
                runs += 1
                if loc[0] == 'F':
                    flagged_runs += 1
                else:
                    bare_runs += 1
            else:
                bare_pair_runs += 1

    for loc in locs:
        for payload, prob in validated_channel(noise, loc).items():
            if prob == 0:
                continue
            if (loc[2], payload) not in allowed[loc[:2]]:
                raise ValueError(f'Tableau fault_list does not cover {loc}, {payload}')
            faults = {loc[:2]: {loc[2]: payload}}
            trial(faults, (*loc, payload))
            if loc[0] == 'B':
                faults[trigger_loc[:2]] = {trigger_loc[2]: 'FLIP'}
                trial(faults, (*loc, payload), require_ft=False)
    for q in DATA:
        for p in 'XYZ':
            trial({}, ('incoming', q, p), incoming=(q, p))
    return {'single_fault_runs': runs, 'flagged_fault_runs': flagged_runs,
            'bare_only_fault_runs': bare_runs, 'incoming_runs': incoming_runs,
            'single_fault_failures': len(failures),
            'cross_check_runs': runs, 'bare_triggered_pair_cross_check_runs': bare_pair_runs,
            'mismatches': len(mismatches), 'failure_details': failures,
            'mismatch_details': mismatches}


def run_preset(params, shots, seed):
    noise = IonNoise(params)
    checks = check_channel_faults(noise)
    print(f'{params.name}: checks {json.dumps(checks)}', flush=True)
    if checks['single_fault_failures'] or checks['mismatches']:
        return {'parameters': asdict(params), 'checks': checks, 'status': 'FAILED'}
    counts = count_faults(noise)
    if counts['single_fault_failures']:
        raise RuntimeError('Fault counting disagrees with single-fault checks; stop')
    print(f'{params.name}: P1={counts["P1"]:.8g}, P2={counts["P2"]:.8g}, '
          f'survival-corrected={counts["P2_survival_corrected"]:.8g}; '
          f'Monte Carlo {shots} shots...', flush=True)
    fails = monte_carlo(noise, shots, random.Random(seed))
    pl = fails / shots
    err = math.sqrt(max(pl, 1 / shots) * (1 - pl) / shots)
    comparisons = {}
    for key in ('P2', 'P2_survival_corrected', 'P2_survival_mean_approx'):
        comparisons[key] = {'estimate_over_MC': counts[key] / pl if pl else None,
                            'estimate_minus_MC_in_stderr': (counts[key] - pl) / err if err else None}
    mc = {'seed': seed, 'shots': shots, 'fails': fails, 'p_L': pl, 'stderr': err}
    print(f'{params.name}: MC={pl:.8g} +/- {err:.3g}; {json.dumps(comparisons)}', flush=True)
    return {'parameters': asdict(params), 'checks': checks, 'counts': counts,
            'monte_carlo': mc, 'comparison': comparisons, 'status': 'PASS'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preset', choices=['ALL', *PRESETS], default='ALL')
    parser.add_argument('--shots', type=int, default=3_000_000)
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'results/stage1_ion_noise.json')
    args = parser.parse_args()
    if args.shots <= 0:
        parser.error('--shots must be positive')
    selected = PRESETS.values() if args.preset == 'ALL' else [PRESETS[args.preset]]
    report = {'stage': 1, 'presets': {}}
    for params in selected:
        result = run_preset(params, args.shots, args.seed)
        report['presets'][params.name] = result
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
        if result['status'] != 'PASS':
            raise SystemExit('Fault-tolerance/cross-check failure: report written; stop for review.')
    print(f'Wrote {args.output}', flush=True)


if __name__ == '__main__':
    main()
