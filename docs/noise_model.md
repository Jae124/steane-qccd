# Noise model

## Stage 0: channel interface (no physical change)

A noise model implements `channel(loc) -> {payload: probability}`, where `loc`
is `(stage, circuit_index, op_index, op)`. Entries are nonidentity fault
outcomes at one location; the remaining probability is no fault. Entries must
be finite and nonnegative, with total at most one (a tolerance of 1e-12 permits
floating-point summation). Unsupported payloads raise an error. An empty channel
means no noise. Explicit zero-probability entries are allowed and remain in
exhaustive enumeration.

The supported payloads and their timing remain those of `steane_ft.fault_list`:
Paulis immediately after resets/CNOTs, and reported-bit `FLIP` at measurements.
Reset Paulis that leave the prepared state unchanged are permitted; they simply
have no physical effect. Both simulators already support all these payloads.

`UniformNoise(p)` uses X after R, Z after RX, FLIP at M/MX, each with probability
p; each nonidentity two-qubit Pauli after CX has probability p/15.
`NoiseModel(p)` and its existing optional reset/measurement rate arguments are
retained. `prob(loc)` remains available for compatibility; the counting and
sampling routines consume `channel(loc)`.

For independent location channels, counting reports

- `P1 = sum_a P(a) fail(a)` over all single faults;
- `P2 = sum_{loc_a < loc_b} sum_{a,b} P(a) P(b) fail(a,b)`;
- `P2_by_logical_type` (X, Z, any), `P2_by_type_pair`, and
  `P2_by_stage_pair`.

X and Z contributions overlap when both logical components fail; `any` counts
each failure once. P2 is the requested weighted pair sum, without no-fault
survival factors. With P1 = 0 it is the leading term when all channel rates
scale together. At finite rates it omits higher-order effects and is not an
exact logical error probability. Faults assigned to circuits skipped by the
adaptive protocol have no effect; the protocol itself is unchanged.

The legacy A coefficients and their breakdowns are returned for the built-in
uniform-rate models, including `count_faults()` (equivalent to counting with
`UniformNoise(1)`). They are normalized independently of p so p=0 is supported.
For uniform noise, P2 = A['any'] p^2 up to floating-point summation error.
For other channel implementations, only probability contributions are reported.

Monte Carlo first samples whether a fault occurs using the channel total, then
samples its payload conditionally. Equal-weight payloads use the original
`rng.choice` path to retain the legacy seeded random sequence. Unequal weights
use inverse cumulative sampling. Latent faults are sampled for every location;
the simulator executes only the circuits selected by the existing protocol.

Assumptions: independent noise across locations; at most one payload per
location; ideal logical state preparation and final decoding, as in the
original simulation. No ion parameters, idle noise, or transport have been
introduced at this checkpoint. The source of the uniform model is the original
repository implementation, rather than a hardware calibration.
