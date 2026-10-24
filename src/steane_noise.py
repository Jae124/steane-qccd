"""
Layer 2 of the thesis stack: logical error rate of one Steane EC round.

Builds on steane_ft.py (layer 1: the flagged EC gadget, verified FT).

What this file does
-------------------
1. A Pauli-frame simulator. Every measurement in a noiseless EC round is
   deterministic (syndrome 0, flag 0), so we only need to track how injected
   Pauli errors propagate. Exact for Clifford circuits + Pauli noise, and
   ~100x faster than a tableau simulator. Cross-checked against the tableau
   simulator from steane_ft.py on every single fault.
2. A per-location noise model. Every fault location (each reset, CNOT and
   measurement, in both the flagged round and the bare round) gets its own
   probability. Today: one global p. Later (layer 3): probabilities derived
   from the QCCD schedule (idle time, shuttling) are dropped in here.
3. Malignant-pair counting: enumerate every pair of faults, count the ones that
   cause a logical error -> exact leading-order coefficient A in p_L ~ A p^2.
4. Monte Carlo at higher p to validate the A p^2 approximation.
5. Plot: p_L vs p for the Steane EC round, against a bare physical qubit (p_L = p).
"""

import itertools
import json
import math
import random
import time

import numpy as np

from steane_ft import (ROUND, SUPPORTS, SUPPORTS_1IDX, stabilizer_circuit,
                       flagged_decode, standard_decode, syndrome_value,
                       fault_list, ec_round, prepare_logical, apply_pauli,
                       output_is_correctable)
import stim

DATA_MASK = (1 << 7) - 1
SUPPORT_MASKS = [sum(1 << q for q in s) for s in SUPPORTS]

# Pre-build the 12 circuits (6 flagged, 6 bare)
CIRCUITS = {
    ('F', i): stabilizer_circuit(kind, SUPPORTS[k], flagged=True)
    for i, (kind, k) in enumerate(ROUND)
}
CIRCUITS.update({
    ('B', i): stabilizer_circuit(kind, SUPPORTS[k], flagged=False)
    for i, (kind, k) in enumerate(ROUND)
})


# ---------------------------------------------------------------------------
# Pauli-frame simulation of one adaptive EC round
# ---------------------------------------------------------------------------

def _apply_payload(fx, fz, qubits, payload):
    for q, p in zip(qubits, payload):
        bit = 1 << q
        if p in 'XY':
            fx ^= bit
        if p in 'ZY':
            fz ^= bit
    return fx, fz


def run_ops_frame(fx, fz, ops, faults):
    """faults: dict op_index -> payload ('XZ', 'X', 'FLIP', ...)."""
    results = {}
    for i, op in enumerate(ops):
        kind = op[0]
        if kind in ('R', 'RX'):
            q = op[1]
            fx &= ~(1 << q)
            fz &= ~(1 << q)
            if i in faults:
                fx, fz = _apply_payload(fx, fz, [q], faults[i])
        elif kind == 'CX':
            c, t = op[1], op[2]
            if (fx >> c) & 1:
                fx ^= 1 << t
            if (fz >> t) & 1:
                fz ^= 1 << c
            if i in faults:
                fx, fz = _apply_payload(fx, fz, [c, t], faults[i])
        elif kind == 'IDLE':
            if i in faults:
                fx, fz = _apply_payload(fx, fz, [op[1]], faults[i])
        else:  # M (Z basis: sees X frame) or MX (X basis: sees Z frame)
            q = op[1]
            bit = (fx >> q) & 1 if kind == 'M' else (fz >> q) & 1
            if faults.get(i) == 'FLIP':
                bit ^= 1
            results[op[2]] = bit
    return fx, fz, results


def ec_round_frame(fx, fz, faults, *, circuits=None, log=None):
    """
    faults: dict (stage, circuit_index) -> {op_index: payload}
    stage 'F' = flagged round, 'B' = bare round (only runs if triggered).
    Returns the data frame after correction.
    """
    circuits = CIRCUITS if circuits is None else circuits
    if log is not None:
        log.update(trigger=None, trigger_index=None)
    trigger = None
    for i, (kind, k) in enumerate(ROUND):
        fx, fz, r = run_ops_frame(fx, fz, circuits[('F', i)], faults.get(('F', i), {}))
        if r['flag'] or r['synd']:
            trigger = (kind, k, r['flag'])
            if log is not None:
                log.update(trigger=trigger, trigger_index=i, trigger_bits=r)
            break
    if trigger is None:
        return fx, fz

    bits = {'Z': [], 'X': []}
    for i, (kind, k) in enumerate(ROUND):
        fx, fz, r = run_ops_frame(fx, fz, circuits[('B', i)], faults.get(('B', i), {}))
        bits[kind].append(r['synd'])
    s_x, s_z = syndrome_value(bits['Z']), syndrome_value(bits['X'])

    if log is not None:
        log.update(s_x=s_x, s_z=s_z)
    kind, k, flag = trigger
    if flag and kind == 'Z':
        z_fix, x_fix = flagged_decode(s_z, SUPPORTS[k]), standard_decode(s_x)
    elif flag and kind == 'X':
        z_fix, x_fix = standard_decode(s_z), flagged_decode(s_x, SUPPORTS[k])
    else:
        z_fix, x_fix = standard_decode(s_z), standard_decode(s_x)
    if log is not None:
        log.update(x_fix=x_fix, z_fix=z_fix)
    for q in x_fix:
        fx ^= 1 << q
    for q in z_fix:
        fz ^= 1 << q
    return fx, fz


def _parity(x):
    return bin(x).count('1') & 1


def logical_failure(fx, fz):
    """
    Ideal decoder on the residual data error, then: did a logical X or Z survive?
    After ideal correction the residual is a stabilizer (even weight) or a
    logical (odd weight), so overall parity decides.
    Returns (x_fail, z_fail).
    """
    dx, dz = fx & DATA_MASK, fz & DATA_MASK
    s_x = syndrome_value([_parity(dx & m) for m in SUPPORT_MASKS])
    s_z = syndrome_value([_parity(dz & m) for m in SUPPORT_MASKS])
    for q in standard_decode(s_x):
        dx ^= 1 << q
    for q in standard_decode(s_z):
        dz ^= 1 << q
    return _parity(dx), _parity(dz)


# ---------------------------------------------------------------------------
# Fault locations and the per-location noise model
# ---------------------------------------------------------------------------

def physical_variants(op):
    """Fault variants at a location under circuit-level depolarizing noise."""
    if op[0] == 'R':
        return ['X']                      # prepared |1> instead of |0>
    if op[0] == 'RX':
        return ['Z']                      # prepared |-> instead of |+>
    if op[0] == 'IDLE':
        return ['Z']
    if op[0] == 'CX':
        return [a + b for a, b in itertools.product('IXYZ', repeat=2) if a + b != 'II']
    return ['FLIP']                       # measurement reports the wrong bit


def all_locations(circuits=None):
    """Every fault location in the gadget: (stage, circuit_index, op_index, op)."""
    circuits = CIRCUITS if circuits is None else circuits
    locs = []
    for stage in ('F', 'B'):
        for i in range(len(ROUND)):
            for j, op in enumerate(circuits[(stage, i)]):
                locs.append((stage, i, j, op))
    return locs


def loc_type(loc):
    op = loc[3][0]
    if op == 'IDLE':
        return 'idle'
    return '2q' if op == 'CX' else ('reset' if op in ('R', 'RX') else 'meas')


class NoiseModel:
    """Legacy operation-rate model, now exposing explicit Pauli channels."""

    def __init__(self, p_2q, p_reset=None, p_meas=None):
        self.p_2q = p_2q
        self.p_reset = p_2q if p_reset is None else p_reset
        self.p_meas = p_2q if p_meas is None else p_meas

    def prob(self, loc):
        op = loc[3]
        if op[0] == 'CX':
            return self.p_2q
        if op[0] in ('R', 'RX'):
            return self.p_reset
        return self.p_meas

    def channel(self, loc):
        variants = physical_variants(loc[3])
        return {v: self.prob(loc) / len(variants) for v in variants}


class UniformNoise(NoiseModel):
    """Each location fails with probability p, uniformly over its variants."""

    def __init__(self, p):
        super().__init__(p)


def validated_channel(noise, loc):
    """Check the nonidentity distribution before counting or sampling it."""
    channel = dict(noise.channel(loc))
    allowed = {payload for _, payload in fault_list([loc[3]])}
    if any(payload not in allowed for payload in channel):
        raise ValueError(f"Unsupported payload at {loc}: {channel}")
    if any(not np.isfinite(p) or p < 0 for p in channel.values()):
        raise ValueError(f"Invalid channel probabilities at {loc}: {channel}")
    if sum(channel.values()) > 1 + 1e-12:
        raise ValueError(f"Channel probability exceeds one at {loc}: {channel}")
    return channel


# ---------------------------------------------------------------------------
# 1. Cross-check frame simulator vs tableau simulator on all single faults
# ---------------------------------------------------------------------------

def cross_check():
    mismatches, n = 0, 0
    for idx in range(len(ROUND)):
        ops = CIRCUITS[('F', idx)]
        for fault in fault_list(ops):
            fx, fz = ec_round_frame(0, 0, {('F', idx): {fault[0]: fault[1]}})
            xf, zf = logical_failure(fx, fz)
            for basis, frame_fail in (('0', xf), ('+', zf)):
                sim = stim.TableauSimulator()
                prepare_logical(sim, basis)
                ec_round(sim, fault_at=(idx, fault))
                tableau_fail = not output_is_correctable(sim, basis)
                n += 1
                mismatches += int(bool(frame_fail) != tableau_fail)
    print(f"Cross-check frame vs tableau simulator: {n} single-fault runs, "
          f"{mismatches} mismatches")
    assert mismatches == 0
    return {"runs": n, "mismatches": mismatches}


# ---------------------------------------------------------------------------
# 2. Single faults and malignant pairs
# ---------------------------------------------------------------------------

def _count_faults_independent(noise=None):
    """Enumerate singles and distinct-location pairs using explicit channels.

    P1/P2 are weighted failure sums, not a finite-rate expansion including
    no-fault factors. A is retained for uniform models and the no-argument API.
    """
    noise = UniformNoise(1.0) if noise is None else noise
    circuits = getattr(noise, 'circuits', CIRCUITS)
    locs = all_locations(circuits)
    channels = [validated_channel(noise, loc) for loc in locs]
    variants = [list(channel) for channel in channels]
    uniform = (type(noise) in (NoiseModel, UniformNoise)
               and noise.p_2q == noise.p_reset == noise.p_meas)
    P1 = 0.0
    P2 = {'X': 0.0, 'Z': 0.0, 'any': 0.0}
    p2_stage, p2_type = {}, {}
    rates = [min(1.0, sum(ch.values())) for ch in channels]
    corrected = 0.0
    corrected_stage, corrected_type = {}, {}

    # First order: any single fault that fails would break fault tolerance
    single_fail = 0
    for loc, vs, channel in zip(locs, variants, channels):
        for v in vs:
            fx, fz = ec_round_frame(0, 0, {(loc[0], loc[1]): {loc[2]: v}}, circuits=circuits)
            if any(logical_failure(fx, fz)):
                single_fail += 1
                P1 += channel[v]

    # Second order: weight each pair by (1/n_variants_i)(1/n_variants_j) so that
    # p_L ~ A p^2 with every location failing with probability p
    A = {'X': 0.0, 'Z': 0.0, 'any': 0.0}
    malignant_pairs = 0
    per_stage = {}
    per_type = {}
    for a, b in itertools.combinations(range(len(locs)), 2):
        la, lb = locs[a], locs[b]
        if not variants[a] or not variants[b]:
            continue
        w = 1.0 / (len(variants[a]) * len(variants[b]))
        survival = math.prod(1 - rate for idx, rate in enumerate(rates)
                             if idx not in (a, b))
        for va in variants[a]:
            for vb in variants[b]:
                faults = {}
                faults.setdefault((la[0], la[1]), {})[la[2]] = va
                faults.setdefault((lb[0], lb[1]), {})[lb[2]] = vb
                fx, fz = ec_round_frame(0, 0, faults, circuits=circuits)
                xf, zf = logical_failure(fx, fz)
                pw = channels[a][va] * channels[b][vb]
                if xf or zf:
                    malignant_pairs += 1
                    A['any'] += w
                    key = la[0] + lb[0]
                    per_stage[key] = per_stage.get(key, 0.0) + w
                    tkey = '*'.join(sorted((loc_type(la), loc_type(lb))))
                    per_type[tkey] = per_type.get(tkey, 0.0) + w
                    P2['any'] += pw
                    corrected += pw * survival
                    corrected_stage[key] = corrected_stage.get(key, 0.0) + pw * survival
                    corrected_type[tkey] = corrected_type.get(tkey, 0.0) + pw * survival
                    p2_stage[key] = p2_stage.get(key, 0.0) + pw
                    p2_type[tkey] = p2_type.get(tkey, 0.0) + pw
                A['X'] += w * xf
                A['Z'] += w * zf
                P2['X'] += pw * xf
                P2['Z'] += pw * zf

    n_loc_types = {}
    for l in locs:
        n_loc_types[l[3][0]] = n_loc_types.get(l[3][0], 0) + 1
    return {
        'n_locations': len(locs),
        'location_counts': n_loc_types,
        'single_fault_failures': single_fail,
        'malignant_pairs': malignant_pairs,
        **({'A': A, 'A_by_stage_pair': per_stage,
            'A_by_type_pair': per_type} if uniform else {}),
        'P1': P1,
        'P2': P2['any'],
        'P2_survival_corrected': corrected,
        'P2_survival_by_stage_pair': corrected_stage,
        'P2_survival_by_type_pair': corrected_type,
        'mean_location_fault_probability': sum(rates) / len(rates),
        'P2_survival_mean_approx': P2['any'] * (1 - sum(rates) / len(rates)) ** (len(locs) - 2),
        'P2_by_logical_type': P2,
        'P2_by_stage_pair': p2_stage,
        'P2_by_type_pair': p2_type,
    }


# ---------------------------------------------------------------------------
# 3. Monte Carlo
# ---------------------------------------------------------------------------

TRIGGER_TYPES = ('reset', 'measurement', '2q', 'idle')


def trigger_contributors(faults, circuits, log):
    """Categories contributing to the nonzero bit(s) at the actual trigger.

    Propagation is linear over GF(2). Replay each category through the executed
    flagged prefix without taking intermediate branches; do not credit faults
    whose contributions cancel or only affect a zero bit at the trigger.
    """
    groups = {}
    stop = log['trigger_index']
    for key, hits in faults.items():
        if key[0] != 'F' or key[1] > stop:
            continue
        for idx, payload in hits.items():
            kind = loc_type((*key, idx, circuits[key][idx]))
            kind = 'measurement' if kind == 'meas' else kind
            groups.setdefault(kind, {}).setdefault(key, {})[idx] = payload
    contributors = []
    for kind, subset in groups.items():
        fx = fz = 0
        for idx in range(stop + 1):
            fx, fz, bits = run_ops_frame(fx, fz, circuits[('F', idx)],
                                        subset.get(('F', idx), {}))
        if any(log['trigger_bits'][tag] and bits[tag] for tag in ('synd', 'flag')):
            contributors.append(kind)
    return contributors


def _monte_carlo_independent(noise, shots, rng, *, return_details=False):
    circuits = getattr(noise, 'circuits', CIRCUITS)
    locs = all_locations(circuits)
    channels = [validated_channel(noise, loc) for loc in locs]
    probs = [min(1.0, sum(ch.values())) for ch in channels]
    variants = [list(ch) for ch in channels]
    equal_weights = [len(set(ch.values())) <= 1 for ch in channels]
    cumulative = [np.cumsum(list(ch.values())).tolist() for ch in channels]
    fails = triggers = bare_fails = 0
    causes = dict.fromkeys((*TRIGGER_TYPES, 'mixed'), 0)
    contributed = dict.fromkeys(TRIGGER_TYPES, 0)
    fractional = dict.fromkeys(TRIGGER_TYPES, 0.0)
    trigger_circuits = [0] * len(ROUND)
    for _ in range(shots):
        faults = {}
        for l, p, vs, equal, edges in zip(locs, probs, variants, equal_weights, cumulative):
            # Preserve the legacy RNG sequence for uniform conditional channels.
            if rng.random() < p:
                if equal:
                    payload = rng.choice(vs)
                else:
                    draw = rng.random() * p
                    payload = next(v for v, edge in zip(vs, edges) if draw < edge)
                faults.setdefault((l[0], l[1]), {})[l[2]] = payload
        if not faults:
            continue
        log = {} if return_details else None
        fx, fz = ec_round_frame(0, 0, faults, circuits=circuits, log=log)
        failed = any(logical_failure(fx, fz))
        fails += int(failed)
        if return_details and log['trigger'] is not None:
            triggers += 1
            trigger_circuits[log['trigger_index']] += 1
            bare_fails += int(failed)
            contributors = trigger_contributors(faults, circuits, log)
            if not contributors:
                raise RuntimeError('A triggered noiseless-input shot has no signal contributor')
            causes[contributors[0] if len(contributors) == 1 else 'mixed'] += 1
            for kind in contributors:
                contributed[kind] += 1
                fractional[kind] += 1 / len(contributors)
    if not return_details:
        return fails
    rate = triggers / shots
    pl = fails / shots
    return {
        'shots': shots, 'fails': fails, 'p_L': pl,
        'stderr': math.sqrt(max(pl, 1 / shots) * (1 - pl) / shots),
        'triggers': triggers, 'trigger_rate': rate,
        'trigger_circuit_counts': trigger_circuits,
        'no_trigger_shots': shots - triggers,
        'trigger_rate_stderr': math.sqrt(rate * (1 - rate) / shots),
        'trigger_cause_counts': causes,
        'trigger_cause_rates': {k: v / shots for k, v in causes.items()},
        'trigger_contributing_type_rates': {k: v / shots for k, v in contributed.items()},
        'trigger_fractional_type_rates': {k: v / shots for k, v in fractional.items()},
        'failures_through_bare': bare_fails,
        'failure_share_through_bare': bare_fails / fails if fails else None,
        'failure_share_through_bare_stderr': (math.sqrt(
            (bare_fails / fails) * (1 - bare_fails / fails) / fails) if fails else None),
    }


def count_faults(noise=None):
    if getattr(noise, 'dynamic_channels', False):
        from ion_transport import conditional_count_faults
        return conditional_count_faults(noise)
    return _count_faults_independent(noise)


def monte_carlo(noise, shots, rng, *, return_details=False):
    if getattr(noise, 'dynamic_channels', False):
        from ion_transport import conditional_monte_carlo
        return conditional_monte_carlo(noise, shots, rng, return_details=return_details)
    return _monte_carlo_independent(noise, shots, rng, return_details=return_details)


if __name__ == '__main__':
    cross_check()

    t0 = time.time()
    counts = count_faults()
    print(f"\nFault counting ({time.time() - t0:.0f} s):")
    print(json.dumps(counts, indent=2))

    rng = random.Random(2026)
    mc = []
    for p in [5e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2, 2e-2, 3e-2, 5e-2, 7e-2, 1e-1]:
        shots = {5e-4: 3_000_000, 1e-3: 1_500_000, 2e-3: 500_000}.get(
            p, 200_000 if p < 1e-2 else 50_000)
        t0 = time.time()
        metrics = monte_carlo(NoiseModel(p), shots, rng, return_details=True)
        f = metrics['fails']
        pl = f / shots
        err = np.sqrt(max(pl, 1 / shots) * (1 - pl) / shots)
        pair_term = counts['A']['any'] * p * p
        mc.append({'p': p, **metrics, 'P2': pair_term,
                   'P2_survival_corrected': pair_term * (1-p)**(counts['n_locations']-2)})
        print(f"MC p={p:.0e}: p_L={pl:.3e} +- {err:.1e}  "
              f"(A p^2 = {counts['A']['any'] * p * p:.3e})  [{time.time() - t0:.0f} s]")

        print(f"  trigger rate={metrics['trigger_rate']:.6g}; "
              f"trigger causes={metrics['trigger_cause_rates']}; "
              f"failure share through bare={metrics['failure_share_through_bare']}")

    with open('../results/steane_ec_results.json', 'w') as fh:
        json.dump({'primary_estimate': 'Monte Carlo', 'counts': counts, 'monte_carlo': mc}, fh, indent=2)
    print("\nWrote results/steane_ec_results.json")
