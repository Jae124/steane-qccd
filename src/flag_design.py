"""Opt-in flag check order and timing policies. Approved instrument is unmodified."""

from functools import partial
import math
from types import FunctionType

import steane_ft as ft
import steane_noise as frame
from ion_noise import IonNoise, IdleIonNoise, check_channel_faults
from ion_timing import idle_probability, insert_idles, operation_duration
from ion_transport import (TransportPlan, conditional_count_faults,
                           conditional_monte_carlo, insert_transport_idles)


def validate_round_order(round_order=None):
    if round_order is None:
        return tuple(ft.ROUND)
    try:
        order = tuple(tuple(entry) for entry in round_order)
        valid = (len(order) == len(ft.ROUND) and set(order) == set(ft.ROUND)
                 and all(len(entry) == 2 and type(entry[1]) is int for entry in order))
    except (TypeError, IndexError):
        valid = False
    if not valid:
        raise ValueError('round_order must be a permutation of the six (kind, support_index) entries of ROUND')
    return order


def order_from_kinds(kinds):
    next_row = {'Z': 0, 'X': 0}
    order = []
    for kind in kinds:
        if kind not in next_row:
            raise ValueError('Expected a six-character Z/X order')
        order.append((kind, next_row[kind]))
        next_row[kind] += 1
    return validate_round_order(order)


CHECK_ORDERS = {name: order_from_kinds(name) for name in ('ZZZXXX', 'XXXZZZ', 'ZXZXZX', 'XZXZXZ')}


def TWO_QUBIT_ONLY(p):
    """CX depolarization only; reset/readout errors and idle noise are zero."""
    return frame.NoiseModel(p, p_reset=0, p_meas=0)


def _decode(trigger, sx, sz, use_flags=True):
    kind, k, flag = trigger
    if use_flags and flag and kind == 'Z':
        return ft.standard_decode(sx), ft.flagged_decode(sz, ft.SUPPORTS[k])
    if use_flags and flag and kind == 'X':
        return ft.flagged_decode(sx, ft.SUPPORTS[k]), ft.standard_decode(sz)
    return ft.standard_decode(sx), ft.standard_decode(sz)


def ec_round_frame(fx, fz, faults, *, circuits=None, log=None, round_order=None):
    order = validate_round_order(round_order)
    if order == tuple(ft.ROUND):
        return frame.ec_round_frame(fx, fz, faults, circuits=circuits, log=log)
    circuits = base_circuits(order) if circuits is None else circuits
    if log is not None:
        log.update(trigger=None, trigger_index=None)
    trigger = None
    for i, (kind, k) in enumerate(order):
        fx, fz, r = frame.run_ops_frame(fx, fz, circuits[('F', i)], faults.get(('F', i), {}))
        if r['flag'] or r['synd']:
            trigger = kind, k, r['flag']
            if log is not None:
                log.update(trigger=trigger, trigger_index=i, trigger_bits=r)
            break
    if trigger is None:
        return fx, fz
    bits = {'Z': [0]*3, 'X': [0]*3}
    for i, (kind, k) in enumerate(order):
        fx, fz, r = frame.run_ops_frame(fx, fz, circuits[('B', i)], faults.get(('B', i), {}))
        bits[kind][k] = r['synd']
    sx, sz = ft.syndrome_value(bits['Z']), ft.syndrome_value(bits['X'])
    x_fix, z_fix = _decode(trigger, sx, sz)
    if log is not None:
        log.update(s_x=sx, s_z=sz, x_fix=x_fix, z_fix=z_fix)
    for q in x_fix:
        fx ^= 1 << q
    for q in z_fix:
        fz ^= 1 << q
    return fx, fz


def ec_round(sim, fault_at=None, use_flags=True, *, faults=None, circuits=None, round_order=None):
    order = validate_round_order(round_order)
    if order == tuple(ft.ROUND):
        return ft.ec_round(sim, fault_at, use_flags, faults=faults, circuits=circuits)
    if fault_at is not None and faults is not None:
        raise ValueError('Specify fault_at or faults, not both')
    injected = faults or {}
    if fault_at is not None:
        i, (op_idx, payload) = fault_at
        injected = {('F', i): {op_idx: payload}}
    circuits = base_circuits(order) if circuits is None else circuits
    log = {'trigger': None, 'trigger_index': None}
    for i, (kind, k) in enumerate(order):
        r = ft.run_ops(sim, circuits[('F', i)], faults=injected.get(('F', i), {}))
        if r['flag'] or r['synd']:
            log.update(trigger=(kind, k, r['flag']), trigger_index=i, trigger_bits=r)
            break
    if log['trigger'] is None:
        return log
    bits = {'Z': [0]*3, 'X': [0]*3}
    for i, (kind, k) in enumerate(order):
        r = ft.run_ops(sim, circuits[('B', i)], faults=injected.get(('B', i), {}))
        bits[kind][k] = r['synd']
    sx, sz = ft.syndrome_value(bits['Z']), ft.syndrome_value(bits['X'])
    x_fix, z_fix = _decode(log['trigger'], sx, sz, use_flags)
    for q in x_fix:
        ft.apply_pauli(sim, [q], 'X')
    for q in z_fix:
        ft.apply_pauli(sim, [q], 'Z')
    log.update(s_x=sx, s_z=sz, x_fix=x_fix, z_fix=z_fix)
    return log


def base_circuits(round_order=None):
    order = validate_round_order(round_order)
    if order == tuple(ft.ROUND):
        return frame.CIRCUITS
    return {(stage, i): ft.stabilizer_circuit(kind, ft.SUPPORTS[k], flagged=stage == 'F')
            for stage in ('F', 'B') for i, (kind, k) in enumerate(order)}


def _family(op):
    return 'reset' if op[0] in ('R', 'RX') else ('meas' if op[0] in ('M', 'MX') else None)


def layer_timing(ops, params, policy='SEQ', events=()):
    """Greedy disjoint adjacent MR layers; transport is a serial layer barrier."""
    if policy not in ('SEQ', 'PAR_MR'):
        raise ValueError('timing_policy must be SEQ or PAR_MR')
    grouped = {}
    for idx, qubit, kind in events:
        grouped.setdefault(idx, []).append((qubit, kind))
    clock = 0.0
    starts, ends = {}, {}
    i = 0
    readout = 0.0
    while i < len(ops):
        clock = math.fsum([clock, *(getattr(params, kind+'_duration') for _, kind in grouped.get(i, []))])
        group = [i]
        if policy == 'PAR_MR' and _family(ops[i]) is not None:
            seen = {ops[i][1]}
            j = i + 1
            while (j < len(ops) and j not in grouped and _family(ops[j]) == _family(ops[i])
                   and ops[j][1] not in seen):
                seen.add(ops[j][1]); group.append(j); j += 1
        spans = []
        for idx in group:
            duration = operation_duration(ops[idx], params)
            starts[idx], ends[idx] = clock, math.fsum([clock, duration])
            if _family(ops[idx]) == 'meas' and params.measurement_duration:
                readout_start = clock + (params.rotation_duration if ops[idx][0] == 'MX' else 0)
                spans.append((readout_start, ends[idx]))
        if spans:
            spans.sort(); lo, hi = spans[0]
            for a, b in spans[1:]:
                if a <= hi: hi = max(hi, b)
                else: readout += hi-lo; lo, hi = a, b
            readout += hi-lo
        clock = max(ends[idx] for idx in group)
        i += len(group)
    clock = math.fsum([clock, *(getattr(params, kind+'_duration') for _, kind in grouped.get(len(ops), []))])
    return starts, ends, clock, readout


def parallel_idles(ops, params, events=()):
    starts, ends, duration, _ = layer_timing(ops, params, 'PAR_MR', events)
    last = {q: 0.0 for q in ft.DATA}
    out, mapping = [], {}
    for idx, op in enumerate(ops):
        touched = op[1:3] if op[0] == 'CX' else (op[1],)
        for q in touched:
            if q in last:
                gap = starts[idx]-last[q]
                if gap > 1e-15:
                    out.append(('IDLE', q, gap))
        mapping[idx] = len(out)
        out.append(op)
        for q in touched:
            last[q] = ends[idx]
            if op[0] in ('M', 'MX'):
                last.pop(q, None)
    for q in ft.DATA:
        gap = duration-last[q]
        if gap > 1e-15:
            out.append(('IDLE', q, gap))
    return out, mapping


def _bind(kernel, **overrides):
    """Private dependency injection: reuse frozen kernel bytecode, no global mutation."""
    namespace = dict(kernel.__globals__)
    namespace.update(overrides)
    result = FunctionType(kernel.__code__, namespace, kernel.__name__, kernel.__defaults__, kernel.__closure__)
    result.__kwdefaults__ = kernel.__kwdefaults__
    return result


class FlagDesign:
    """Experimental view of a noise model with explicit order/policy and unchanged decoder."""
    def __init__(self, noise, round_order=None, *, timing_policy='SEQ', transport_events=None, heating_penalty=None):
        self.round_order = validate_round_order(round_order)
        if timing_policy not in ('SEQ', 'PAR_MR'):
            raise ValueError('timing_policy must be SEQ or PAR_MR')
        self.noise = noise
        self.params = getattr(noise, 'params', None)
        self.timed = isinstance(noise, IdleIonNoise)
        if not self.timed and timing_policy != 'SEQ':
            raise ValueError('PAR_MR requires a timed preset; H2/uniform/two-qubit-only remain timing-free')
        self.timing_policy = timing_policy if self.timed else 'timing-free'
        self.base_circuits = base_circuits(self.round_order)
        inherited = getattr(noise, 'transport_events', ())
        self.transport_events = tuple(inherited if transport_events is None else transport_events)
        penalties = getattr(noise, 'heating_penalty', None) if heating_penalty is None else heating_penalty
        self.transport_plan = None
        self.dynamic_channels = False
        if self.transport_events:
            if not self.timed:
                raise ValueError('Transport requires sourced timing')
            inserter = parallel_idles if timing_policy == 'PAR_MR' else None
            self.transport_plan = TransportPlan(self.params, self.transport_events, penalties,
                                                circuits=self.base_circuits, idle_inserter=inserter)
            self.circuits = self.transport_plan.circuits
            self.dynamic_channels = self.transport_plan.dynamic_channels
        elif not self.timed or math.isinf(self.params.T2):
            self.circuits = self.base_circuits
        elif timing_policy == 'SEQ':
            self.circuits = (noise.circuits if self.round_order == tuple(ft.ROUND) and not inherited else
                             {key: insert_idles(ops, self.params) for key, ops in self.base_circuits.items()})
        else:
            self.circuits = {key: parallel_idles(ops, self.params)[0] for key, ops in self.base_circuits.items()}
        from ion_transport import validate_penalties
        same_penalties = not self.transport_events or validate_penalties(penalties) == getattr(noise, 'heating_penalty', {})
        self.legacy = (self.round_order == tuple(ft.ROUND) and timing_policy == 'SEQ'
                       and self.transport_events == tuple(inherited) and same_penalties)
        self._frame = partial(ec_round_frame, round_order=self.round_order)
        self._tableau = partial(ec_round, round_order=self.round_order)

    def channel(self, loc):
        if self.transport_plan and loc[3][0] == 'CX':
            return {v: self.transport_plan.cx_rate(loc)/15 for v in frame.physical_variants(loc[3])}
        if self.timed and loc[3][0] == 'IDLE':
            p = idle_probability(loc[3][2], self.params.T2)
            return {'Z': p} if p else {}
        # The operation channels are independent of their circuit's position.
        if self.timed:
            return IonNoise.channel(self.noise, loc)
        return self.noise.channel(loc)

    def channel_for_trigger(self, loc, trigger_index):
        if self.transport_plan and loc[3][0] == 'CX':
            return {v: self.transport_plan.cx_rate(loc, trigger_index)/15 for v in frame.physical_variants(loc[3])}
        return self.channel(loc)

    def count_faults(self):
        if self.legacy:
            return frame.count_faults(self.noise)
        kernel = conditional_count_faults if self.dynamic_channels else frame._count_faults_independent
        return _bind(kernel, ec_round_frame=self._frame)(self)

    def monte_carlo(self, shots, rng, *, return_details=True):
        if self.legacy:
            return frame.monte_carlo(self.noise, shots, rng, return_details=return_details)
        kernel = conditional_monte_carlo if self.dynamic_channels else frame._monte_carlo_independent
        return _bind(kernel, ec_round_frame=self._frame)(self, shots, rng, return_details=return_details)

    def check_faults(self):
        if self.legacy:
            return check_channel_faults(self.noise)
        return _bind(check_channel_faults, ec_round_frame=self._frame, ec_round=self._tableau)(self)

    def exposure(self, mc):
        if not self.timed:
            return None
        from ion_timing import idle_exposure
        if self.legacy:
            plan = getattr(self.noise, 'transport_plan', None)
            return (idle_exposure(self.params, mc, circuits=self.noise.circuits, extra_durations=plan.extra_durations)
                    if plan else idle_exposure(self.params, mc))
        events = self.transport_plan.grouped if self.transport_plan else {}
        times = {key: layer_timing(ops, self.params, self.timing_policy, events.get(key, ()))[2:]
                 for key, ops in self.base_circuits.items()}
        paths = {}
        for stop in [None, *range(6)]:
            keys = ([('F', i) for i in range(6)] if stop is None else
                    [('F', i) for i in range(stop+1)] + [('B', i) for i in range(6)])
            count = mc['no_trigger_shots'] if stop is None else mc['trigger_circuit_counts'][stop]
            paths['no_trigger' if stop is None else f'trigger_F{stop}'] = {
                'count': count, 'probability': count/mc['shots'],
                'round_duration_seconds': math.fsum(times[key][0] for key in keys),
                'readout_idle_seconds_per_data_qubit': math.fsum(times[key][1] for key in keys),
                'data': {str(q): {'seconds': math.fsum(o[2] for key in keys for o in self.circuits[key]
                                                     if o[0]=='IDLE' and o[1]==q),
                                  'summed_Z_probability': math.fsum(idle_probability(o[2],self.params.T2)
                                         for key in keys for o in self.circuits[key] if o[0]=='IDLE' and o[1]==q)}
                         for q in ft.DATA}}
        average = {}
        for q in ft.DATA:
            average[str(q)] = {}
            for metric in ('seconds', 'summed_Z_probability'):
                mean = math.fsum(p['probability']*p['data'][str(q)][metric] for p in paths.values())
                var = math.fsum(p['probability']*(p['data'][str(q)][metric]-mean)**2 for p in paths.values())
                average[str(q)][metric] = mean
                average[str(q)][metric+'_stderr'] = math.sqrt(var/mc['shots'])
        return {'definition': 'Actual executed paths under the selected order and timing policy.',
                'shots': mc['shots'], 'paths': paths, 'monte_carlo_weighted_data': average,
                'mean_readout_idle_seconds_per_data_qubit': math.fsum(
                    p['probability']*p['readout_idle_seconds_per_data_qubit'] for p in paths.values())}
