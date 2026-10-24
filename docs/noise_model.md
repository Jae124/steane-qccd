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

## Stage 1: operation-specific trapped-ion channels

All physical inputs live in `src/ion_params.py`, in immutable `IonParams`
presets. Stage 1 uses probabilities only; durations and T2 do not cause noise
yet. No circuit gate order, decoder rule, or adaptive stopping rule changed.
The tableau simulator now also accepts the frame simulator's injection mapping,
including multiple locations and faults in the bare follow-up round.

### Inputs and provenance

`JONES_MURALI_2026` keeps the requested preset name. Its source is Jones and
Murali, [arXiv:2510.23519v1, Table 1 and Section 5.1](https://arxiv.org/html/2510.23519v1#S5),
posted October 2025. It supplies T2 = 2.2 s; MS/rotation/measurement/reset
40/5/400/50 microseconds; shuttle/split/merge/junction 5/80/80/100 microseconds;
measurement probability 0.001 and reset probability 0.005. The two-qubit 0.001
and one-qubit 0.00001 probabilities are **study placeholders**, not evaluations
of the paper's heating-dependent fidelity model. A logical CX is treated as
one effective noisy two-qubit gate at the MS error rate; a native gate
synthesis and its additional rotations are outside Stage 1.

`H2_LIKE` uses the initialization/one-qubit/two-qubit/measurement values
0.00004/0.000029/0.00128/0.001 from LeBlond et al.,
[arXiv:2508.14227v1, Table 1](https://arxiv.org/html/2508.14227v1#S3.T1).
That paper attributes the gate and initialization values to the H2-1E emulator
specification; its measurement value is a compromise between biased measurement
rates. We use these numbers as requested, directly as total nonidentity Pauli
probabilities. Their correspondence to average or process infidelity is
**unverified**; no fidelity conversion is silently applied. Timing and T2
are unspecified (`None`) for this preset because Stage 1 does not use them.
A later idle-noise stage must explicitly choose and document those inputs;
we have not assigned H2 the Jones–Murali timing by implication.

The optional `fidelity_to_depolarizing(F_avg, n_qubits)` helper uses
`p = (d+1)/d * (1-F_avg)` for `d = 2**n_qubits`: factors 3/2 and 5/4 for one
and two qubits. Here p is total probability of a nonidentity Pauli, not the
maximally-mixed replacement parameter. This follows from
`F_avg = (d*F_ent + 1)/(d+1)` and `F_ent = 1-p` for a Pauli channel. Inputs
that would give p outside [0, 1] are rejected.

### Exact local channel composition

Let `r = p_reset`, `g = p_1q`, and `m = p_meas`.

- CX: each of the 15 nonidentity Paulis has probability `p_2q/15`.
- R: X with probability r.
- RX: prepare |+> with a Z error of probability r (the rotated reset bit flip),
  then apply independent depolarizing rotation noise. Multiplying the Paulis
  modulo phase gives `P(X) = P(Y) = g/3`,
  `P(Z) = (1-r)*g/3 + r*(1-g)`. The omitted identity probability is
  `(1-r)*(1-g) + r*g/3`. This includes cancellation exactly. An X outcome
  on |+> has no physical effect, but remains an explicit Pauli outcome.
- M: FLIP with probability m.
- MX: marginalize the depolarizing rotation immediately before X measurement.
  X commutes with that measurement; Y and Z anticommute, giving rotation-flip
  probability `b = 2*g/3`. Combining independently with the readout flip gives
  `P(FLIP) = m + b - 2*m*b` by XOR, including double-flip cancellation.

The MX representation is exact for the joint reported measurement record and
the remaining data/ancilla state. The measured ion is unused until reset (or
round end), so its postmeasurement state can be discarded. No new op or fault
payload is needed: both simulators already support the effective FLIP. The
unit tests compare explicit premeasurement Paulis with this reduction on every
check circuit containing MX, both logical bases, and an incoming data error.
This reduction must be revisited if a future circuit reuses a measured ion
without resetting it. RX/MX combine multiple hardware error mechanisms into
one effective location; the pair expansion and N below refer to these effective
locations. Rotation contributions are grouped under reset or measurement
according to the operation they implement; no standalone 1q gate is present.

### Fault coverage and diagnostics

The ion channel support is R: X; RX: X/Y/Z; CX: 15 Paulis; M/MX: FLIP.
The existing tableau `fault_list` already covers that entire support. Per
preset, testing both logical bases gives 1,152 flagged fault runs and 756
bare-only runs (1,908 total), plus 42 incoming-error trials. Bare-only faults
are skipped by the noiseless flagged round. To actually exercise the bare
injections, another 756 cross-check trials force the first flagged syndrome
measurement to flip, then inject each bare fault. Those are two-fault agreement
tests, not single-fault tolerance claims. Cross-checking compares six stabilizer
expectations and the relevant logical observable before ideal output decoding,
as well as the logical failure decision. Any fault-tolerance failure or
simulator mismatch is saved with its location, trigger, syndromes, correction,
and residual frame, and stops the CLI before estimation.

### Survival correction (required for this and every later stage)

For each location let `p_loc = sum(channel(loc).values())`. The counter now
also reports the exact requested full-location survival-weighted pair sum:

```
P2_survival_corrected = sum_{i<j,a,b} P_i(a) P_j(b) fail(a,b)
                       * product_{k != i,j} (1-p_k)
```

The product is computed independently for each pair, including all 96 potential
locations in flagged and bare circuits, whether executed or skipped. Direct
products avoid division by zero when a location has probability one. Corrected
breakdowns by type pair and stage pair sum to this corrected total.
`P2_survival_mean_approx = P2*(1-p_bar)**(N-2)` is also included for comparison.
For uniform channels both corrections equal `P2*(1-p)**94` to rounding error.

This exact product gives the logical failure probability from exactly two
nonidentity assignments in the **full latent location pool**. It is a subset
of Monte Carlo failures, not a complete finite-rate estimate: Monte Carlo also
includes three or more assignments, including assignments to skipped circuits
that have no operational effect. Requiring survival on skipped locations
therefore suppresses this term even when those extra faults are irrelevant.
Combining RX/MX mechanisms also means it counts effective-location faults,
not primitive hardware faults. The mean-rate approximation has neither the
exact pair weighting nor a general lower-bound guarantee. Raw P2 can over- or
underestimate full Monte Carlo because it also omits survival and higher orders.
Reports give both estimates, their ratios to Monte Carlo, and differences in
Monte Carlo standard errors; disagreement does not imply a tolerance failure.
