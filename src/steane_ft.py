"""
Steane [[7,1,3]] code: flagged syndrome extraction (Chao-Reichardt style)
+ exhaustive single-fault injection test.

Qubit layout (0-indexed in code, 1-indexed in comments/tables):
    data qubits 0..6  (= qubits 1..7 of the code)
    7 = syndrome ancilla 'a'
    8 = flag ancilla     'f'

Protocol for one EC round (adaptive):
  1. Flagged round: measure the 6 stabilizers one by one with flagged circuits.
     Stop early as soon as a flag fires or a syndrome bit is 1.
  2. If triggered: one bare (unflagged) round of all 6 stabilizers.
  3. Decode: flagged lookup table for the hooked error type if a flag fired,
     standard (syndrome = qubit address) otherwise.

The test: inject every possible single fault into step 1, run the protocol,
then check the output data has at most a correctable (weight <= 1) error.
That is the definition of a fault-tolerant distance-3 EC gadget.
"""

import itertools
import stim

DATA = list(range(7))
A, F = 7, 8
N_QUBITS = 9

# Supports of the three stabilizer rows (1-indexed), from the Hamming matrix H
SUPPORTS_1IDX = [(4, 5, 6, 7), (2, 3, 6, 7), (1, 3, 5, 7)]
SUPPORTS = [tuple(q - 1 for q in s) for s in SUPPORTS_1IDX]


# ---------------------------------------------------------------------------
# Circuits, written as lists of operations so we can inject faults anywhere
# ---------------------------------------------------------------------------

def stabilizer_circuit(kind, support, flagged=True):
    """
    kind = 'Z': measure Z_a Z_b Z_c Z_d  (detects X errors)
         ancilla a in |0>, CNOT data -> a, measure a in Z
         flag f in |+>, CNOT f -> a after 1st and before last data CNOT, measure f in X
    kind = 'X': measure X_a X_b X_c X_d  (detects Z errors)
         Hadamard-dual: ancilla in |+>, CNOT a -> data, measure a in X
         flag f in |0>, CNOT a -> f at the same two places, measure f in Z
    Returns a list of ops: ('R', q) | ('RX', q) | ('CX', c, t) | ('M', q, tag) | ('MX', q, tag)
    """
    q1, q2, q3, q4 = support
    ops = []
    if kind == 'Z':
        ops.append(('R', A))
        if flagged:
            ops.append(('RX', F))
        data_cx = lambda q: ('CX', q, A)
        flag_cx = ('CX', F, A)
    else:
        ops.append(('RX', A))
        if flagged:
            ops.append(('R', F))
        data_cx = lambda q: ('CX', A, q)
        flag_cx = ('CX', A, F)

    ops.append(data_cx(q1))
    if flagged:
        ops.append(flag_cx)
    ops.append(data_cx(q2))
    ops.append(data_cx(q3))
    if flagged:
        ops.append(flag_cx)
    ops.append(data_cx(q4))

    if kind == 'Z':
        ops.append(('M', A, 'synd'))
        if flagged:
            ops.append(('MX', F, 'flag'))
    else:
        ops.append(('MX', A, 'synd'))
        if flagged:
            ops.append(('M', F, 'flag'))
    return ops


# Order of the 6 stabilizer measurements in a round: Z-checks then X-checks
ROUND = [('Z', k) for k in range(3)] + [('X', k) for k in range(3)]


def ops_to_stim(ops):
    """Convert our op list into a stim.Circuit (for drawing / export)."""
    c = stim.Circuit()
    for op in ops:
        if op[0] in ('R', 'RX', 'M', 'MX'):
            c.append(op[0], [op[1]])
        elif op[0] == 'CX':
            c.append('CX', [op[1], op[2]])
        c.append('TICK')
    return c


# ---------------------------------------------------------------------------
# Execution with an optional injected fault
# ---------------------------------------------------------------------------

def apply_pauli(sim, qubits, paulis):
    ps = stim.PauliString(N_QUBITS)
    for q, p in zip(qubits, paulis):
        ps[q] = p
    sim.do(ps)


def run_ops(sim, ops, fault=None):
    """
    Execute ops on a TableauSimulator.
    fault = (op_index, payload): after op op_index, apply a Pauli (payload is a
    string like 'XZ' over the op's qubits) or, for a measurement, flip the
    result (payload == 'FLIP').
    Returns dict tag -> measured bit.
    """
    results = {}
    for i, op in enumerate(ops):
        kind = op[0]
        hit = fault is not None and fault[0] == i
        if kind == 'R':
            sim.reset(op[1])
        elif kind == 'RX':
            sim.reset_x(op[1])
        elif kind == 'CX':
            sim.cnot(op[1], op[2])
        elif kind in ('M', 'MX'):
            q = op[1]
            if kind == 'MX':
                sim.h(q)
            bit = sim.measure(q)
            if kind == 'MX':
                sim.h(q)
            if hit and fault[1] == 'FLIP':
                bit = not bit
            results[op[2]] = int(bit)
            continue
        if hit:
            qubits = [op[1]] if kind in ('R', 'RX') else [op[1], op[2]]
            apply_pauli(sim, qubits, fault[1])
    return results


def fault_list(ops):
    """Every single fault for a circuit: Paulis after resets / CNOTs, flipped measurements."""
    faults = []
    for i, op in enumerate(ops):
        if op[0] in ('R', 'RX'):
            faults += [(i, p) for p in 'XYZ']
        elif op[0] == 'CX':
            faults += [(i, p1 + p2) for p1, p2 in itertools.product('IXYZ', repeat=2)
                       if (p1, p2) != ('I', 'I')]
        elif op[0] in ('M', 'MX'):
            faults.append((i, 'FLIP'))
    return faults


# ---------------------------------------------------------------------------
# Decoders
# ---------------------------------------------------------------------------

def standard_decode(s):
    """Syndrome s (0..7) -> list of data qubits (0-indexed) to flip."""
    return [] if s == 0 else [s - 1]


def flagged_decode(s, support):
    """
    Flagged lookup table for a circuit with CNOT order (a, b, c, d):
      0 -> none, a -> a, d -> d, c^d -> {c, d}; anything else -> standard
    (Syndromes and addresses here are 1-indexed qubit numbers.)
    """
    a, b, c, d = (q + 1 for q in support)
    table = {0: [], a: [a], d: [d], c ^ d: [c, d]}
    if s in table:
        return [q - 1 for q in table[s]]
    return standard_decode(s)


def syndrome_value(bits):
    """Three bits (row r1, r2, r3) -> binary number = qubit address."""
    return 4 * bits[0] + 2 * bits[1] + bits[2]


# ---------------------------------------------------------------------------
# One full adaptive EC round
# ---------------------------------------------------------------------------

def ec_round(sim, fault_at=None, use_flags=True):
    """
    fault_at = (circuit_index in ROUND, fault) or None.
    Returns a log dict describing what happened.
    """
    log = {'trigger': None}
    trigger = None
    for idx, (kind, k) in enumerate(ROUND):
        ops = stabilizer_circuit(kind, SUPPORTS[k], flagged=True)
        fault = fault_at[1] if (fault_at and fault_at[0] == idx) else None
        r = run_ops(sim, ops, fault)
        if r['flag'] or r['synd']:
            trigger = (kind, k, r['flag'])
            break
    log['trigger'] = trigger
    if trigger is None:
        return log

    # Bare round (fault-free by the single-fault assumption)
    bits = {'Z': [], 'X': []}
    for kind, k in ROUND:
        r = run_ops(sim, stabilizer_circuit(kind, SUPPORTS[k], flagged=False))
        bits[kind].append(r['synd'])
    s_x = syndrome_value(bits['Z'])  # Z-checks see X errors
    s_z = syndrome_value(bits['X'])  # X-checks see Z errors

    kind, k, flag = trigger
    # A flagged Z-check circuit can hook Z errors; a flagged X-check circuit hooks X errors
    if use_flags and flag and kind == 'Z':
        z_fix, x_fix = flagged_decode(s_z, SUPPORTS[k]), standard_decode(s_x)
    elif use_flags and flag and kind == 'X':
        z_fix, x_fix = standard_decode(s_z), flagged_decode(s_x, SUPPORTS[k])
    else:
        z_fix, x_fix = standard_decode(s_z), standard_decode(s_x)

    for q in x_fix:
        apply_pauli(sim, [q], 'X')
    for q in z_fix:
        apply_pauli(sim, [q], 'Z')
    log.update(s_x=s_x, s_z=s_z, x_fix=x_fix, z_fix=z_fix)
    return log


# ---------------------------------------------------------------------------
# State prep and checks
# ---------------------------------------------------------------------------

def prepare_logical(sim, basis):
    """Ideal (noiseless) encoder for |0_L>; |+_L> = H^7 |0_L>."""
    for q in (0, 1, 3):                  # pivots: qubits 1, 2, 4
        sim.h(q)
    for pivot, targets in ((3, (4, 5, 6)), (1, (2, 5, 6)), (0, (2, 4, 6))):
        for t in targets:
            sim.cnot(pivot, t)
    if basis == '+':
        for q in DATA:
            sim.h(q)


def pauli_on(support, p):
    ps = stim.PauliString(N_QUBITS)
    for q in support:
        ps[q] = p
    return ps


def output_is_correctable(sim, basis):
    """
    Ideal decoder on the output: read the syndromes noiselessly, apply the
    standard correction, check the logical observable. Passes iff the residual
    error had weight <= 1 (up to stabilizers).
    """
    s_x = syndrome_value([int(sim.peek_observable_expectation(pauli_on(s, 'Z')) == -1)
                          for s in SUPPORTS])
    s_z = syndrome_value([int(sim.peek_observable_expectation(pauli_on(s, 'X')) == -1)
                          for s in SUPPORTS])
    for q in standard_decode(s_x):
        apply_pauli(sim, [q], 'X')
    for q in standard_decode(s_z):
        apply_pauli(sim, [q], 'Z')
    logical = pauli_on(DATA, 'Z' if basis == '0' else 'X')
    return sim.peek_observable_expectation(logical) == 1


def check_encoder():
    for basis in ('0', '+'):
        sim = stim.TableauSimulator()
        prepare_logical(sim, basis)
        for s in SUPPORTS:
            assert sim.peek_observable_expectation(pauli_on(s, 'Z')) == 1
            assert sim.peek_observable_expectation(pauli_on(s, 'X')) == 1
        obs = pauli_on(DATA, 'Z' if basis == '0' else 'X')
        assert sim.peek_observable_expectation(obs) == 1
    print("Encoder check: |0_L> and |+_L> are +1 eigenstates of all 6 stabilizers and the logical. OK")


# ---------------------------------------------------------------------------
# The single-fault injection test
# ---------------------------------------------------------------------------

def single_fault_test(use_flags=True, verbose_failures=5):
    # Faults inside the 6 flagged circuits
    locations = []
    for idx, (kind, k) in enumerate(ROUND):
        ops = stabilizer_circuit(kind, SUPPORTS[k], flagged=True)
        locations += [(idx, f) for f in fault_list(ops)]
    # Plus a single error already on the data when the round starts
    incoming = [(q, p) for q in DATA for p in 'XYZ']

    total, failures = 0, []
    for basis in ('0', '+'):            # |0_L> catches X_L, |+_L> catches Z_L
        for loc in locations:
            sim = stim.TableauSimulator()
            prepare_logical(sim, basis)
            log = ec_round(sim, fault_at=loc, use_flags=use_flags)
            total += 1
            if not output_is_correctable(sim, basis):
                failures.append((basis, loc, log))
        for q, p in incoming:
            sim = stim.TableauSimulator()
            prepare_logical(sim, basis)
            apply_pauli(sim, [q], p)
            log = ec_round(sim, use_flags=use_flags)
            total += 1
            if not output_is_correctable(sim, basis):
                failures.append((basis, ('incoming', q + 1, p), log))

    label = "WITH flagged decoding" if use_flags else "standard decoder only (flags ignored)"
    print(f"\nSingle-fault test, {label}:")
    print(f"  {total} single-fault runs, {len(failures)} logical failures")
    for basis, loc, log in failures[:verbose_failures]:
        if loc[0] == 'incoming':
            where = f"incoming {loc[2]} on qubit {loc[1]}"
        else:
            kind, k = ROUND[loc[0]]
            op = stabilizer_circuit(kind, SUPPORTS[k])[loc[1][0]]
            where = (f"{kind}-check on {SUPPORTS_1IDX[k]}, after op #{loc[1][0]} "
                     f"{op}, fault {loc[1][1]}")
        print(f"   - |{basis}_L>: {where} -> trigger {log['trigger']}, "
              f"s_x={log.get('s_x')}, s_z={log.get('s_z')}, "
              f"X fix on {[q + 1 for q in log.get('x_fix', [])]}, "
              f"Z fix on {[q + 1 for q in log.get('z_fix', [])]}")
    return len(failures)


if __name__ == '__main__':
    check_encoder()
    n_fail = single_fault_test(use_flags=True)
    single_fault_test(use_flags=False)
    print("\nRESULT:", "flagged EC round is fault-tolerant (0 failures)" if n_fail == 0
          else "NOT fault-tolerant, see failures above")

    # Export one flagged circuit and one bare circuit as Stim circuits + diagrams
    flagged = ops_to_stim(stabilizer_circuit('Z', SUPPORTS[0], flagged=True))
    with open('../results/flagged_Z4567.stim', 'w') as fh:
        fh.write(str(flagged))
    with open('../figures/flagged_Z4567.svg', 'w') as fh:
        fh.write(str(flagged.diagram('timeline-svg')))
    full = stim.Circuit()
    for kind, k in ROUND:
        full += ops_to_stim(stabilizer_circuit(kind, SUPPORTS[k], flagged=True))
    with open('../results/flagged_round.stim', 'w') as fh:
        fh.write(str(full))
    print("\nWrote results/flagged_*.stim and figures/flagged_Z4567.svg")
