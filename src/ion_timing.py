"""Sequential chain timing, with exact dephasing intervals cut at branch boundaries."""

import math

from steane_ft import DATA
from steane_noise import CIRCUITS


def idle_probability(duration, T2):
    if not math.isfinite(duration) or duration < 0:
        raise ValueError('Idle duration must be finite and nonnegative')
    if T2 is None or math.isnan(T2) or T2 <= 0:
        raise ValueError('A positive sourced T2 is required')
    return -math.expm1(-duration / T2) / 2


def operation_duration(op, params):
    fields = {'CX': ('ms_duration',), 'R': ('reset_duration',),
              'RX': ('reset_duration', 'rotation_duration'),
              'M': ('measurement_duration',),
              'MX': ('rotation_duration', 'measurement_duration')}
    if op[0] not in fields:
        raise ValueError(f'Unsupported timed operation {op}')
    values = [getattr(params, key) for key in fields[op[0]]]
    if any(value is None or not math.isfinite(value) or value < 0 for value in values):
        raise ValueError(f'Sourced nonnegative durations are required for {op[0]}')
    return sum(values)


def insert_idles(ops, params):
    """One Z channel per uninterrupted local wait, closed at circuit boundaries.

    Data are live throughout this interval, including before first/after last
    participation. Each ancilla is live only from reset completion to the start
    of its terminal measurement. No physical time is added by an IDLE marker.
    """
    # Compute intervals first so summation roundoff cannot create tiny idle gaps.
    duration = [operation_duration(op, params) for op in ops]
    last = {q: -1 for q in DATA}
    live = set(DATA)
    out = []
    for idx, op in enumerate(ops):
        touched = op[1:3] if op[0] == 'CX' else (op[1],)
        for q in touched:
            if q in live:
                wait = math.fsum(duration[last[q] + 1:idx])
                if wait > 0:
                    out.append(('IDLE', q, wait))
        out.append(op)
        for q in touched:
            last[q] = idx
            if op[0] in ('R', 'RX'):
                live.add(q)
            elif op[0] in ('M', 'MX'):
                live.discard(q)
    for q in DATA:
        wait = math.fsum(duration[last[q] + 1:])
        if wait > 0:
            out.append(('IDLE', q, wait))
    return out


def timed_circuits(params):
    # Validate even the infinite-T2 path; timing-free presets stay unsupported.
    idle_probability(0.0, params.T2)
    for ops in CIRCUITS.values():
        for op in ops:
            operation_duration(op, params)
    if math.isinf(params.T2):
        # Exact regression, including RNG stream and location numbering.
        return CIRCUITS
    return {key: insert_idles(ops, params) for key, ops in CIRCUITS.items()}
