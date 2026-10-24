"""Optional event hooks, including heating carried along actual adaptive paths.

Event indices refer to original stabilizer_circuit ops, before idle insertion.
Nothing here generates or asserts a hardware transport schedule.
"""

import itertools
import math

import numpy as np

from ion_params import DEFAULT_HEATING_PENALTY
from ion_timing import operation_duration
from steane_ft import DATA, N_QUBITS, ROUND
from steane_noise import (CIRCUITS, TRIGGER_TYPES, all_locations, ec_round_frame,
                          loc_type, logical_failure, run_ops_frame,
                          trigger_contributors, validated_channel)


def validate_penalties(penalties=None):
    result = dict(DEFAULT_HEATING_PENALTY)
    if penalties is not None:
        from collections.abc import Mapping
        if not isinstance(penalties, Mapping):
            raise ValueError('heating_penalty must be a mapping')
        if set(penalties) - set(result):
            raise ValueError('Unknown heating penalty kind')
        result.update(penalties)
    try:
        valid = all(math.isfinite(value) and value >= 0 for value in result.values())
    except (TypeError, ValueError):
        valid = False
    if not valid:
        raise ValueError('Heating penalties must be finite and nonnegative')
    return result


def validate_events(params, events, *, circuits=None):
    circuits = CIRCUITS if circuits is None else circuits
    result = []
    for event in events:
        if not isinstance(event, (list, tuple)) or len(event) != 5:
            raise ValueError('Transport event must be (stage, circuit, before_op, qubit, kind)')
        stage, circuit, before, qubit, kind = event
        if stage not in ('F', 'B') or type(circuit) is not int or not 0 <= circuit < len(ROUND):
            raise ValueError(f'Invalid event circuit: {event}')
        if type(before) is not int or not 0 <= before <= len(circuits[(stage, circuit)]):
            raise ValueError(f'Invalid original before_op index: {event}')
        if type(qubit) is not int or not 0 <= qubit < N_QUBITS:
            raise ValueError(f'Invalid moved qubit: {event}')
        if kind not in DEFAULT_HEATING_PENALTY:
            raise ValueError(f'Invalid transport kind: {event}')
        duration = getattr(params, kind + '_duration')
        if duration is None or not math.isfinite(duration) or duration < 0:
            raise ValueError(f'A sourced finite duration is required for {kind}')
        result.append((stage, circuit, before, qubit, kind))
    return tuple(result)


def insert_transport_idles(ops, params, events):
    """Accumulate all waiting, including transport on every currently live ion."""
    grouped = {}
    for before, qubit, kind in events:
        grouped.setdefault(before, []).append((qubit, kind))
    waiting = {q: [] for q in DATA}
    live = set(DATA)
    out, mapping = [], {}
    for idx in range(len(ops) + 1):
        for _, kind in grouped.get(idx, []):
            for q in live:
                waiting[q].append(getattr(params, kind + '_duration'))
        if idx == len(ops):
            break
        op = ops[idx]
        touched = op[1:3] if op[0] == 'CX' else (op[1],)
        for q in touched:
            if q in live:
                gap = math.fsum(waiting[q])
                if gap:
                    out.append(('IDLE', q, gap))
                waiting[q] = []
        mapping[idx] = len(out)
        out.append(op)
        duration = operation_duration(op, params)
        for q in live - set(touched):
            waiting[q].append(duration)
        for q in touched:
            if op[0] in ('R', 'RX'):
                live.add(q)
                waiting[q] = []
            elif op[0] in ('M', 'MX'):
                live.discard(q)
                waiting.pop(q, None)
    for q in DATA:
        gap = math.fsum(waiting[q])
        if gap:
            out.append(('IDLE', q, gap))
    return out, mapping


class TransportPlan:
    def __init__(self, params, transport_events, heating_penalty=None, *, circuits=None, idle_inserter=None):
        self.params = params
        base_circuits = CIRCUITS if circuits is None else circuits
        idle_inserter = insert_transport_idles if idle_inserter is None else idle_inserter
        self.events = validate_events(params, transport_events, circuits=base_circuits)
        self.penalties = validate_penalties(heating_penalty)
        self.grouped = {}
        for stage, circuit, before, qubit, kind in self.events:
            self.grouped.setdefault((stage, circuit), []).append((before, qubit, kind))
        self.circuits, self.mapping, self.extra_durations = {}, {}, {}
        for key, ops in base_circuits.items():
            events = self.grouped.get(key, [])
            self.circuits[key], self.mapping[key] = idle_inserter(ops, params, events)
            self.extra_durations[key] = math.fsum(getattr(params, kind + '_duration')
                                                 for _, _, kind in events)
        self.rates = {}
        for stop in [None, *range(len(ROUND))]:
            # None supplies the virtual bare table after a full flagged prefix;
            # these bare faults remain unexecuted in no-trigger shots.
            keys = [('F', i) for i in range(len(ROUND) if stop is None else stop + 1)]
            keys += [('B', i) for i in range(len(ROUND))]
            pending = [0.0] * N_QUBITS
            table = {}
            for key in keys:
                by_index = {}
                for before, qubit, kind in self.grouped.get(key, []):
                    by_index.setdefault(before, []).append((qubit, kind))
                for idx in range(len(base_circuits[key]) + 1):
                    for qubit, kind in by_index.get(idx, []):
                        pending[qubit] += self.penalties[kind]
                    if idx == len(base_circuits[key]):
                        continue
                    op = base_circuits[key][idx]
                    if op[0] == 'CX':
                        a, b = op[1:]
                        rate = params.p_2q + pending[a] + pending[b]
                        if rate > 1:
                            raise ValueError('Heating makes a CX probability exceed one')
                        table[(*key, self.mapping[key][idx])] = rate
                        pending[a] = pending[b] = 0.0
            self.rates[stop] = table
        self.dynamic_channels = any(
            any(table[key] != self.rates[None][key] for key in table if key[0] == 'B')
            for table in self.rates.values())

    def cx_rate(self, loc, trigger_index=None):
        # Flagged prefixes have identical causal histories in all branches.
        table = self.rates[None] if loc[0] == 'F' else self.rates[trigger_index]
        return table[loc[:3]]


def synthetic_ancilla_events():
    """Hand-written-style in/out markers per check; a test fixture, not routing."""
    from steane_ft import A
    result = []
    for stage in ('F', 'B'):
        for circuit in range(len(ROUND)):
            ops = CIRCUITS[(stage, circuit)]
            first = next(i for i, op in enumerate(ops) if op[0] == 'CX')
            readout = next(i for i, op in enumerate(ops)
                           if op[0] in ('M', 'MX') and op[2] == 'synd')
            result += [(stage, circuit, first, A, 'shuttle'),
                       (stage, circuit, readout, A, 'shuttle')]
    return result


class _ContextNoise:
    def __init__(self, noise, trigger):
        self.noise, self.trigger = noise, trigger
    def channel(self, loc):
        return self.noise.channel_for_trigger(loc, self.trigger)


def _context_channels(noise, locs):
    tables = {stop: [validated_channel(_ContextNoise(noise, stop), loc) for loc in locs]
              for stop in [None, *range(len(ROUND))]}
    for i in range(len(locs)):
        support = list(dict.fromkeys(v for table in tables.values() for v in table[i]))
        for table in tables.values():
            for payload in support:
                table[i].setdefault(payload, 0.0)
    return tables


def conditional_count_faults(noise):
    """Pair attribution with B rates conditioned on the executed flagged prefix."""
    locs = all_locations(noise.circuits)
    channels = _context_channels(noise, locs)
    variants = [list(ch) for ch in channels[None]]
    rates = {stop: [min(1.0, sum(ch.values())) for ch in table]
             for stop, table in channels.items()}
    P1 = 0.0
    single = malignant = 0
    P2 = dict.fromkeys(('X', 'Z', 'any'), 0.0)
    stage, types, corrected_stage, corrected_types = {}, {}, {}, {}
    corrected = approximate = 0.0
    for i, loc in enumerate(locs):
        for payload in variants[i]:
            log = {}
            frame = ec_round_frame(0, 0, {loc[:2]: {loc[2]: payload}}, circuits=noise.circuits, log=log)
            if any(logical_failure(*frame)):
                single += 1
                P1 += channels[log['trigger_index']][i][payload]
    for a, b in itertools.combinations(range(len(locs)), 2):
        if not variants[a] or not variants[b]:
            continue
        la, lb = locs[a], locs[b]
        survival = {stop: math.prod(1 - rate for idx, rate in enumerate(rs) if idx not in (a, b))
                    for stop, rs in rates.items()}
        for va in variants[a]:
            for vb in variants[b]:
                faults = {}
                faults.setdefault(la[:2], {})[la[2]] = va
                faults.setdefault(lb[:2], {})[lb[2]] = vb
                log = {}
                frame = ec_round_frame(0, 0, faults, circuits=noise.circuits, log=log)
                xf, zf = logical_failure(*frame)
                if not (xf or zf):
                    continue
                stop = log['trigger_index']
                weight = channels[stop][a][va] * channels[stop][b][vb]
                malignant += 1
                P2['any'] += weight
                P2['X'] += weight * xf
                P2['Z'] += weight * zf
                skey = la[0] + lb[0]
                tkey = '*'.join(sorted((loc_type(la), loc_type(lb))))
                stage[skey] = stage.get(skey, 0.0) + weight
                types[tkey] = types.get(tkey, 0.0) + weight
                value = weight * survival[stop]
                corrected += value
                corrected_stage[skey] = corrected_stage.get(skey, 0.0) + value
                corrected_types[tkey] = corrected_types.get(tkey, 0.0) + value
                mean = sum(rates[stop]) / len(locs)
                approximate += weight * (1-mean)**(len(locs)-2)
    counts = {}
    for loc in locs:
        counts[loc[3][0]] = counts.get(loc[3][0], 0) + 1
    return {'n_locations': len(locs), 'location_counts': counts,
            'single_fault_failures': single, 'malignant_pairs': malignant,
            'P1': P1, 'P2': P2['any'], 'P2_by_logical_type': P2,
            'P2_by_stage_pair': stage, 'P2_by_type_pair': types,
            'P2_survival_corrected': corrected, 'P2_survival_by_stage_pair': corrected_stage,
            'P2_survival_by_type_pair': corrected_types, 'P2_survival_mean_approx': approximate,
            'mean_location_fault_probability': sum(rates[None])/len(locs),
            'mean_location_probability_by_trigger': {str(k): sum(v)/len(v) for k, v in rates.items()},
            'channel_context': 'Bare CX rates and survival factors conditioned on actual flagged stop.'}


def conditional_monte_carlo(noise, shots, rng, *, return_details=False):
    """Sample F, determine the branch, then sample the corresponding B channel."""
    locs = all_locations(noise.circuits)
    contexts = _context_channels(noise, locs)
    sampling = {}
    for stop, channels in contexts.items():
        sampling[stop] = [(loc, min(1.0, sum(ch.values())), list(ch),
                           len(set(ch.values())) <= 1, np.cumsum(list(ch.values())).tolist())
                          for loc, ch in zip(locs, channels)]
    fails = triggers = bare_fails = 0
    causes = dict.fromkeys((*TRIGGER_TYPES, 'mixed'), 0)
    contributed = dict.fromkeys(TRIGGER_TYPES, 0)
    fractional = dict.fromkeys(TRIGGER_TYPES, 0.0)
    trigger_circuits = [0] * len(ROUND)

    def sample(rows, stage, faults):
        for loc, p, variants, equal, edges in rows:
            if loc[0] != stage:
                continue
            if rng.random() < p:
                if equal:
                    payload = rng.choice(variants)
                else:
                    draw = rng.random() * p
                    payload = next(v for v, edge in zip(variants, edges) if draw < edge)
                faults.setdefault(loc[:2], {})[loc[2]] = payload

    for _ in range(shots):
        faults = {}
        sample(sampling[None], 'F', faults)
        stop = None
        if faults:
            fx = fz = 0
            for idx in range(len(ROUND)):
                fx, fz, bits = run_ops_frame(fx, fz, noise.circuits[('F', idx)], faults.get(('F', idx), {}))
                if bits['flag'] or bits['synd']:
                    stop = idx
                    break
        sample(sampling[stop], 'B', faults)
        if not faults:
            continue
        log = {} if return_details else None
        frame = ec_round_frame(0, 0, faults, circuits=noise.circuits, log=log)
        failed = any(logical_failure(*frame))
        fails += int(failed)
        if return_details and log['trigger'] is not None:
            triggers += 1
            trigger_circuits[log['trigger_index']] += 1
            bare_fails += int(failed)
            contributors = trigger_contributors(faults, noise.circuits, log)
            if not contributors:
                raise RuntimeError('A triggered shot has no signal contributor')
            causes[contributors[0] if len(contributors) == 1 else 'mixed'] += 1
            for kind in contributors:
                contributed[kind] += 1
                fractional[kind] += 1 / len(contributors)
    if not return_details:
        return fails
    rate, pl = triggers/shots, fails/shots
    return {'shots': shots, 'fails': fails, 'p_L': pl,
            'stderr': math.sqrt(max(pl, 1/shots)*(1-pl)/shots),
            'triggers': triggers, 'trigger_rate': rate, 'trigger_circuit_counts': trigger_circuits,
            'no_trigger_shots': shots-triggers, 'trigger_rate_stderr': math.sqrt(rate*(1-rate)/shots),
            'trigger_cause_counts': causes, 'trigger_cause_rates': {k:v/shots for k,v in causes.items()},
            'trigger_contributing_type_rates': {k:v/shots for k,v in contributed.items()},
            'trigger_fractional_type_rates': {k:v/shots for k,v in fractional.items()},
            'failures_through_bare': bare_fails,
            'failure_share_through_bare': bare_fails/fails if fails else None,
            'failure_share_through_bare_stderr': (math.sqrt((bare_fails/fails)*(1-bare_fails/fails)/fails)
                                                 if fails else None)}


def run_stage4(*, shots=3_000_000, seed=2026, output=None, transport_events=None,
               heating_penalty=None, synthetic=False):
    """Checkpoint runner: empty default, optional explicitly synthetic fixture."""
    import json
    import hashlib
    from dataclasses import asdict
    from pathlib import Path
    from ion_analysis import fingerprint, save_json
    from ion_noise import IdleIonNoise, run_preset
    from ion_params import JONES_MURALI_2026, SYNTHETIC_HEATING_PENALTIES
    from ion_timing import idle_exposure
    root = Path(__file__).resolve().parents[1]
    output = root / 'results/stage4_transport_results.json' if output is None else Path(output)
    if shots <= 0:
        raise ValueError('shots must be positive')
    if synthetic and (transport_events is not None or heating_penalty is not None):
        raise ValueError('Choose the synthetic fixture or custom event/penalty inputs')
    # Reject invalid custom inputs before spending time on the empty regression.
    if not synthetic:
        IdleIonNoise(JONES_MURALI_2026, transport_events, heating_penalty)
    approved = json.loads((root/'results/ion_noise_results.json').read_text())
    prior = approved['presets']['JONES_MURALI_2026']
    baseline_seed = prior['monte_carlo']['seed']
    report = {'stage': 4, 'primary_estimate': 'Monte Carlo',
              'approved_stage3_sha256': hashlib.sha256((root/'results/ion_noise_results.json').read_bytes()).hexdigest(),
              'schedule_scope': ('Synthetic interface tests only; no hardware schedule.' if synthetic
                                 else 'Explicit event input only; no schedule generated.'),
              'presets': {'H2_LIKE': approved['presets']['H2_LIKE']}, 'synthetic_cases': []}
    report['presets']['H2_LIKE'].update(analysis_stage=4, transport_events=[],
        transport_scope='Stage 1 only; no sourced timing, no transport applied.')
    print('Empty transport regression against approved Stage 3', flush=True)
    baseline = run_preset(JONES_MURALI_2026, shots, baseline_seed, idle=True, transport_events=[])
    report['presets']['JONES_MURALI_2026'] = baseline
    save_json(output, report)
    if baseline['status'] != 'PASS':
        raise RuntimeError('Empty-transport fault check failed; stop for review')
    base_model = IdleIonNoise(JONES_MURALI_2026)
    empty_model = IdleIonNoise(JONES_MURALI_2026, [])
    same_channels = (base_model.circuits == empty_model.circuits and
        all(base_model.channel(l) == empty_model.channel(l) for l in all_locations(base_model.circuits)))
    checks = {'circuits_and_channels_equal': same_channels,
              'counts_equal': baseline['counts'] == prior['counts'],
              'approved_full_MC_compared': shots == prior['monte_carlo']['shots'],
              'MC_equal': (baseline['monte_carlo'] == prior['monte_carlo']
                           if shots == prior['monte_carlo']['shots'] else None),
              'approved_Stage3_cache_compatible': fingerprint() == approved['simulation_fingerprint']}
    report['empty_transport_regression'] = checks
    if not same_channels or not checks['counts_equal'] or checks['MC_equal'] is False or not checks['approved_Stage3_cache_compatible']:
        save_json(output, report)
        raise RuntimeError('Empty transport did not reproduce Stage 3; stop for review')
    baseline.update(stage=4, scope='Stage 3 model with empty transport', transport_events=[],
                    heating_penalty=dict(DEFAULT_HEATING_PENALTY),
                    p_eff=prior['p_eff'], p_eff_P2_secondary=prior['p_eff_P2_secondary'])
    baseline['idle_exposure'] = idle_exposure(JONES_MURALI_2026, baseline['monte_carlo'])
    if synthetic:
        if transport_events is not None or heating_penalty is not None:
            raise ValueError('Choose the synthetic fixture or custom event/penalty inputs')
        events = synthetic_ancilla_events()
        cases = [({'shuttle': value}, seed+4000+i) for i,value in enumerate(SYNTHETIC_HEATING_PENALTIES)]
    else:
        events = list(transport_events or [])
        cases = [(validate_penalties(heating_penalty), seed+4000)] if events else []
    report['event_fixture'] = events
    for penalties, point_seed in cases:
        print(f'Transport hook test: {len(events)} events; penalties={penalties}', flush=True)
        result = run_preset(JONES_MURALI_2026, shots, point_seed, idle=True,
                            transport_events=events, heating_penalty=penalties)
        result.update(stage=4, schedule_scope='Synthetic test' if synthetic else 'User-supplied events',
                      transport_events=events, heating_penalty=validate_penalties(penalties))
        report['synthetic_cases'].append(result)
        save_json(output, report)
        if result['status'] != 'PASS':
            raise RuntimeError('Transport fault-tolerance/cross-check failure; stop for review')
        model = IdleIonNoise(JONES_MURALI_2026, events, penalties)
        result['dynamic_channels'] = model.dynamic_channels
        result['extra_duration_per_potential_circuit'] = [
            {'stage':key[0],'circuit_index':key[1],'seconds':value}
            for key,value in model.transport_plan.extra_durations.items()]
        result['idle_exposure'] = idle_exposure(model.params,result['monte_carlo'],
            circuits=model.circuits,extra_durations=model.transport_plan.extra_durations)
        result['heating_affected_CX'] = [
            {'stage':key[0],'circuit_index':key[1],'timed_op_index':key[2],'p_2q':rate}
            for key,rate in model.transport_plan.rates[None].items() if rate != model.params.p_2q]
        save_json(output, report)
    if synthetic:
        values=[v['counts']['P2'] for v in report['synthetic_cases']]
        report['P2_strictly_increases_with_heating'] = all(a<b for a,b in zip(values,values[1:]))
        if not report['P2_strictly_increases_with_heating']:
            save_json(output, report)
            raise RuntimeError('Synthetic P2 monotonicity failed; stop for review')
    save_json(output, report)
    print(f'Wrote {output}',flush=True)
    return report
