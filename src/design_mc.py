"""Sparse, batched Monte Carlo for the A2 independent-location flag designs.

A Binomial(n,p) count followed by a uniform subset is exactly n independent
Bernoulli location draws. Conditional payloads use the original channel weights.
All nonzero shots run the full adaptive Pauli protocol; zero-fault shots are
known successes. Four component frames retain the original trigger attribution.
This sampler has its own NumPy stream; approved scalar RNG streams are untouched.
"""
import math
import numpy as np

import steane_ft as ft
import steane_noise as scalar

CATEGORIES = ('reset', 'measurement', '2q', 'idle')
LOGICAL = np.array([scalar.logical_failure(i, 0)[0] for i in range(128)], dtype=np.uint8)


def _mask(qubits, payload, letters):
    return sum(1 << q for q, letter in zip(qubits, payload) if letter in letters)


class BatchSampler:
    def __init__(self, design):
        if design.transport_events or design.dynamic_channels:
            raise ValueError('A2 batch sampling supports independent channels and empty transport only')
        self.design = design
        self.locations = scalar.all_locations(design.circuits)
        self.channels = [scalar.validated_channel(design, loc) for loc in self.locations]
        self.by_circuit = {}
        self.compiled = []
        for i, (loc, channel) in enumerate(zip(self.locations, self.channels)):
            op = loc[3]
            qubits = op[1:3] if op[0] == 'CX' else [op[1]]
            variants = list(channel)
            xs = [0] + [_mask(qubits, v, 'XY') if v != 'FLIP' else 0 for v in variants]
            zs = [0] + [_mask(qubits, v, 'ZY') if v != 'FLIP' else 0 for v in variants]
            flips = [0] + [int(v == 'FLIP') for v in variants]
            category = scalar.loc_type(loc)
            category = 'measurement' if category == 'meas' else category
            self.compiled.append((CATEGORIES.index(category), np.array(xs, np.uint16),
                                  np.array(zs, np.uint16), np.array(flips, np.uint8)))
            self.by_circuit.setdefault(loc[:2], []).append(i)
        self.standard = np.array([_mask(ft.standard_decode(s), 'X'*7, 'X') for s in range(8)], np.uint16)
        self.flagged = {k:np.array([sum(1 << q for q in ft.flagged_decode(s, ft.SUPPORTS[k]))
                                   for s in range(8)], np.uint16) for k in range(3)}

    def sample_faults(self, shots, rng):
        """Sparse arrays (shot indices, 1-based channel-payload codes) per location."""
        faults = {}
        for i, channel in enumerate(self.channels):
            probabilities = np.array(list(channel.values()), dtype=float)
            p = float(probabilities.sum())
            count = rng.binomial(shots, min(1.0, p))
            if not count:
                continue
            indices = rng.choice(shots, count, replace=False)
            payloads = rng.choice(len(probabilities), count, p=probabilities/p).astype(np.uint8) + 1
            faults[i] = (indices, payloads)
        return faults

    def simulate(self, shots, faults, *, traces=False):
        """Replay explicit faults; traces allow an exact shot-by-shot scalar oracle."""
        indices = np.unique(np.concatenate([x[0] for x in faults.values()])) if faults else np.array([], int)
        n = len(indices)
        codes = np.zeros((len(self.locations), n), dtype=np.uint8)
        for loc, (hits, payloads) in faults.items():
            codes[loc, np.searchsorted(indices, hits)] = payloads
        fx = np.zeros((4, n), np.uint16); fz = np.zeros_like(fx)
        stop = np.full(n, -1, np.int8); flags = np.zeros(n, np.uint8)
        causes = np.zeros((4, n), bool)
        trigger_synd = np.zeros(n, np.uint8)

        def run(key, active):
            x, z = fx[:, active].copy(), fz[:, active].copy()
            bits = {}
            for loc in self.by_circuit[key]:
                op = self.locations[loc][3]
                kind, q = op[:2]
                category, xs, zs, flips = self.compiled[loc]
                code = codes[loc, active]
                if kind in ('R', 'RX'):
                    x &= np.uint16(65535 ^ (1 << q)); z &= np.uint16(65535 ^ (1 << q))
                elif kind == 'CX':
                    t = op[2]
                    x ^= ((x >> q) & 1) << t
                    z ^= ((z >> t) & 1) << q
                elif kind in ('M', 'MX'):
                    bit = ((x if kind == 'M' else z) >> q) & 1
                    bit[category] ^= flips[code]
                    bits[op[2]] = bit
                    continue
                x[category] ^= xs[code]; z[category] ^= zs[code]
            fx[:, active], fz[:, active] = x, z
            return bits

        for i in range(6):
            active = np.flatnonzero(stop < 0)
            if not len(active): break
            bits = run(('F', i), active)
            synd = np.bitwise_xor.reduce(bits['synd'], axis=0)
            flag = np.bitwise_xor.reduce(bits['flag'], axis=0)
            hit = (synd | flag).astype(bool)
            selected = active[hit]
            stop[selected] = i; flags[selected] = flag[hit]; trigger_synd[selected] = synd[hit]
            causes[:, selected] = ((bits['synd'][:, hit] & synd[hit]) | (bits['flag'][:, hit] & flag[hit])).astype(bool)
        active = np.flatnonzero(stop >= 0)
        sx = np.zeros(n, np.uint8); sz = np.zeros(n, np.uint8)
        for i, (kind, k) in enumerate(self.design.round_order):
            if not len(active): break
            bits = run(('B', i), active)
            synd = np.bitwise_xor.reduce(bits['synd'], axis=0).astype(np.uint8)
            (sx if kind == 'Z' else sz)[active] |= synd << (2-k)
        x = np.bitwise_xor.reduce(fx, axis=0); z = np.bitwise_xor.reduce(fz, axis=0)
        xfix = self.standard[sx]; zfix = self.standard[sz]
        for i, (kind, k) in enumerate(self.design.round_order):
            selected = (stop == i) & (flags != 0)
            if kind == 'Z': zfix[selected] = self.flagged[k][sz[selected]]
            else: xfix[selected] = self.flagged[k][sx[selected]]
        x ^= xfix; z ^= zfix
        xfail, zfail = LOGICAL[x & 127].astype(bool), LOGICAL[z & 127].astype(bool)
        failed = xfail | zfail
        multiplicity = causes.sum(axis=0)
        if np.any((stop >= 0) & (multiplicity == 0)):
            raise RuntimeError('A trigger has no contributing signal type')
        exclusive = {c:int(np.count_nonzero(causes[i] & (multiplicity == 1))) for i,c in enumerate(CATEGORIES)}
        exclusive['mixed'] = int(np.count_nonzero(multiplicity > 1))
        result = {'shots':int(shots), 'fails':int(failed.sum()),
                  'triggers':len(active), 'trigger_circuit_counts':[int(np.count_nonzero(stop == i)) for i in range(6)],
                  'trigger_cause_counts':exclusive,
                  'trigger_contributing_type_counts':{c:int(causes[i].sum()) for i,c in enumerate(CATEGORIES)},
                  'trigger_fractional_type_counts':{c:float(np.sum(causes[i]/np.maximum(multiplicity,1))) for i,c in enumerate(CATEGORIES)},
                  'failures_through_bare':int(np.count_nonzero(failed & (stop >= 0))),
                  'logical_failure_counts':{'X_only':int(np.count_nonzero(xfail & ~zfail)),
                                            'Z_only':int(np.count_nonzero(zfail & ~xfail)),
                                            'both_XZ':int(np.count_nonzero(xfail & zfail))}}
        if traces:
            result['traces'] = {'shot_indices':indices,'fx':x,'fz':z,'stop':stop,'flag':flags,
                                'trigger_synd':trigger_synd,'s_x':sx,'s_z':sz,'x_fix':xfix,'z_fix':zfix,
                                'contributors':causes,'x_fail':xfail,'z_fail':zfail}
        return result

    def sample(self, shots, rng):
        return self.simulate(shots, self.sample_faults(shots, rng))


def merge_counts(total, addition):
    """Integer-count accumulation; no averaging of independent chunk estimates."""
    if not total:
        return addition
    for key, value in addition.items():
        if isinstance(value, dict):
            for label, count in value.items(): total[key][label] += count
        elif isinstance(value, list):
            total[key] = [a+b for a,b in zip(total[key],value)]
        else: total[key] += value
    return total


def metrics(counts):
    n = counts['shots']; failed = counts['fails']; pl = failed/n
    trigger = counts['triggers']/n; bare = counts['failures_through_bare']
    share = bare/failed if failed else None
    return dict(counts, p_L=pl, stderr=math.sqrt(max(pl,1/n)*(1-pl)/n),
                no_trigger_shots=n-counts['triggers'], trigger_rate=trigger,
                trigger_rate_stderr=math.sqrt(trigger*(1-trigger)/n),
                trigger_cause_rates={k:v/n for k,v in counts['trigger_cause_counts'].items()},
                trigger_contributing_type_rates={k:v/n for k,v in counts['trigger_contributing_type_counts'].items()},
                trigger_fractional_type_rates={k:v/n for k,v in counts['trigger_fractional_type_counts'].items()},
                failure_share_through_bare=share,
                failure_share_through_bare_stderr=math.sqrt(share*(1-share)/failed) if failed else None,
                logical_failure_rates={k:v/n for k,v in counts['logical_failure_counts'].items()},
                logical_failure_shares={k:v/failed if failed else None for k,v in counts['logical_failure_counts'].items()})
