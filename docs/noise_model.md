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

## Stage 2: sequential idle dephasing and reporting decisions

Stage 2 applies timing **only to JONES_MURALI_2026**. H2_LIKE remains timing-free,
with results marked **Stage 1 only**, pending a sourced timing preset. Its
unspecified timing fields remain None. No transport schedule or heating model
is introduced.

### Timing and idle placement

`src/ion_timing.py` assumes one chain executing one operation at a time.
Durations come exclusively from `ion_params.py`: CX uses one MS duration;
RX uses reset plus rotation; MX uses rotation plus measurement; R and M use
their own durations. An operation's own participating ions are active for its
whole duration and receive the Stage 1 operation channel. All other live ions
wait. This is an operation-level timing abstraction, without overlapping gates
or pulse-level subdivisions of RX/MX.

Data are live throughout each executed circuit/round. Ancillas become live at
their reset and stop being live at their terminal measurement. No idle noise is
charged before preparation or after measurement. Gates and measurements are
ideal apart from their specified local channels and waiting noise.

Each uninterrupted wait on a qubit within a circuit becomes one
`('IDLE', qubit, duration)` marker immediately before its next operation. At
circuit end, outstanding data waits are closed with trailing markers, including
full-circuit waits on data outside that stabilizer's support. Markers add no
time and apply no noiseless operation. They export as identity gates in Stim;
the durations are retained in the Python ops, not in Stim's TICK instruction.
The only permitted idle payload is Z, in both simulators and `fault_list`.

**Boundary subdivision:** The fixed location representation cuts an untouched
data gap at each circuit boundary, which is a possible adaptive termination
point. The suffix of an executed circuit and the prefix of its successor are
both charged. No part of a skipped circuit is charged in the execution.
Thus a maximal physical gap can have more than one marker across boundaries;
within a circuit, gaps are lumped rather than inserting noise after every op.
For Z dephasing this boundary subdivision is exact, not a timing approximation:

```
p_idle(t) = (1 - exp(-t/T2))/2
p_idle(t1+t2) = p_idle(t1) + p_idle(t2) - 2*p_idle(t1)*p_idle(t2)
```

Z commutes with all intervening operations not touching the ion, so delaying a
wait's fault to its marker has the same action. Tests compare injection at the
beginning and end of every generated gap and verify time conservation for
every data ion in every circuit. This subdivision keeps fixed stage/circuit/op
fault IDs and independent local channels across adaptive branches. The primary
Monte Carlo channel is exact under this abstraction. As with any splitting or
combining of effective locations, raw and survival-corrected P2 attribution
refer to this particular location representation; they are not universal
finite-rate observables.

`T2 = infinity` removes idle markers entirely, reproducing Stage 1's location
indices, all pair-count outputs and seeded sampling exactly. Noiseless runs
with finite-T2 markers preserve trivial syndromes, no flags and the encoded
logical state. Both flagged and bare circuits have waiting locations.

The 400 microsecond measurement wait at T2=2.2 s has p_idle ≈ 9.0901e-5. A whole
gap can contain multiple operations, so the largest probability need not equal
that single-operation example. In this representation the largest gap is a
whole flagged circuit (1150 microseconds) for a data ion outside its support.
The result file records that gap and its probability. The requested T2 sweep
(0.1, 0.5, 2.2, 10 seconds) is a sensitivity study, with values centralized in
`ion_params.py`; it does not assert measured device parameters at every point.

### Monte Carlo and trigger metrics (all stages/presets)

Monte Carlo is the primary p_L estimate from Stage 2 onward. P2, with and
without survival factors, is retained for attribution/comparison only. In
Stage 3, p_eff must be computed from Monte Carlo, with a P2-derived value as
a secondary check; the Stage 3 implementation is described below.

`monte_carlo(..., return_details=True)` reports p_L, its binomial standard
error, trigger rate (bare-round execution probability), its standard error,
and the fraction of failed shots that executed the bare round. The latter is
`failures_through_bare / fails`, with None when there were no failures. Failure
rate samples the union of logical X/Z failure events, as before. Optional
reporting does not consume random numbers or change the sampled faults.

Trigger attribution is **signal provenance at the observed first trigger**,
not a claim of a unique counterfactual cause in every multifault shot. Propagate
each fault category separately through the executed flagged prefix, suppressing
intermediate stopping while evaluating that contribution. Pauli-frame
propagation, resets, and measured-bit flips are linear over GF(2). Credit a
category if its contribution is nonzero in at least one of the syndrome/flag
bits that is actually 1 at the trigger. Ignore bare-round and later flagged
faults. Ignore category contributions that cancel completely or affect only a
zero bit at the actual trigger.

Reports provide three complementary splits:

- `trigger_cause_rates`: reset, measurement, 2q, idle, plus mixed. A shot with
  one contributing category goes to that category; more than one goes to mixed.
  These mutually exclusive rates sum to the trigger rate.
- `trigger_contributing_type_rates`: include each category that contributes;
  these overlap in mixed shots and need not sum to the trigger rate.
- `trigger_fractional_type_rates`: divide mixed shots equally among their
  contributing categories. These four rates sum to the trigger rate. Equal
  sharing is an explicit reporting convention, not uniquely identified causality.

RX/MX mechanisms retain their Stage 1 grouping: an RX rotation fault is in
reset; an MX rotation fault is in measurement. Counts of single-fault trials
and failures explicitly separate flagged locations, bare-only locations
(skipped on noiseless input), and incoming errors. Forced-trigger bare tests
are cross-checks, not added to the single-fault FT headline.

The Stage 2 output also reruns Stage 1 Jones–Murali and timing-free H2-like with
these metrics. `results/stage0_trigger_metrics.json` backfills the same metrics
for the Stage 0 uniform benchmark without changing its exact 542/3M result.
Historical Stage 0/1 result files remain as originally approved.

## Stage 3: effective uniform rate from Monte Carlo

The primary definition is the uniform per-location p for which the **uniform
Monte Carlo** logical failure probability equals the ion Monte Carlo p_L.
`src/ion_analysis.py` builds the calibration, readout sweep, figures and final
JSON; `src/ion_calibration.py` implements interpolation and statistical inference.
`python ion_noise.py --stage 3` (now the default) runs/resumes this analysis.
It adds no new noise channels or transport events.

### Calibration, pseudothreshold and confidence intervals

The uniform grid is centralized in `ion_params.py` and contains 11 p values
from 0.00085 to 0.00185, each with 1,000,000 independent shots. It is the union
of two six-point neighborhoods, sharing p=0.00135, around the timing-free H2
and nominal Jones–Murali targets. Each target's result explicitly lists six
bracketing grid p values and their shot counts. The same grid also brackets
every readout-sensitivity target and the uniform pseudothreshold. Grid and
new target runs use distinct deterministic random seeds, recorded in the JSON.

For neighboring grid points, interpolate linearly in (log p, log p_L). Invert
that piecewise curve to get p_eff; never extrapolate. The uniform
pseudothreshold is the first upward crossing of p_L(p)=p, solving
`log p_L - log p = 0` within the corresponding log-linear segment. This is a
one-gadget pseudothreshold, not a many-round accuracy threshold.

The reported **95% confidence intervals** use 20,000 parametric binomial
bootstrap draws. In each draw, resample each grid failure count at its observed
binomial rate and independently resample each target failure count, then
repeat the log-log interpolation/inversion. Use the same sampled grid for
p_eff and the pseudothreshold, so their covariance is retained in the
p_eff/pseudothreshold ratio interval. p_eff/p_2q uses the preset's fixed p_2q.
The point estimate uses the observed MC means, not an average of bootstrap
roots. Percentile bounds are 2.5% and 97.5% of valid draws.

An inverse curve must be monotone. The observed grid is checked for reversals;
when a bootstrap draw has a small adjacent reversal, weighted pool-adjacent-
violators regression on log rates restores nondecreasing order, with weights
`N*p_L/(1-p_L)` (inverse delta-method log-rate variance). Exact plateau targets
use the rightmost p endpoint. The result records the adjustment fraction and
whether the point curve was adjusted. If more than 0.5% of bootstrap draws
lack a bracket, the analysis stops and requires a wider grid rather than
silently clipping the interval. The unadjusted point grid is strictly
increasing in the saved experiment.

Intervals include target and calibration **sampling** uncertainty conditional
on the log-linear interpolation model. They do not include interpolation bias,
hardware-rate uncertainty, or uncertainty in placeholder gate parameters.
Finite grid resolution is visible in the calibration figure. Unit tests cover
analytic power-law inversion, exact uniform round trips, monotone pooling,
bootstrap reproducibility, and the separate contributions of target and grid
sampling errors.

The secondary value remains `p_eff_P2 = sqrt(P2/A_uniform)`, with the approved
uniform leading coefficient. Its ratio to the primary p_eff and membership in
the primary 95% interval are reported. Survival-corrected pair sums, raw pair
sums and comparisons with MC remain in every target and grid record. Neither
pair quantity substitutes for primary MC p_L or the MC p_eff.

### Readout-duration sensitivity

Jones–Murali measurement duration is set to 0, 100, 200 or 400 microseconds;
every other parameter, including measurement error probability and the 5
microsecond MX basis rotation, stays nominal. Each point has 3,000,000 shots,
its own seed, pair attribution, trigger metrics, separate fault counts and a
calibrated MC p_eff with a propagated interval. The 400 microsecond sample is
a fresh nominal Stage 3 run, rather than reusing the Stage 2 target sample.
Stage 2 files and their T2 sweep remain archived reference data. H2 primary
MC is reused from the approved 3M-shot timing-free Stage 2 report and remains
marked Stage 1 only; no H2 timing or idle exposure is invented.

Zero-duration readout is a counterfactual timing intervention, not a claim of
instantaneous hardware measurement. It removes readout waiting for every
other live ion, including flag/ancilla waits as well as all data waits. Thus
the logical-rate difference isolates readout-duration effects under this
model, but is not an intervention confined exclusively to data qubits.
The report includes the difference between 400 and 0 microseconds, its
sampling uncertainty, and its ratio to the total nominal idle increment over
the approved Stage 1 no-idle reference. That ratio uses a separate binomial
bootstrap of the three target experiments; it is a ratio of logical-rate
increments, not a microscopic assignment of every failure to a mechanism.

### Nominal idle exposure on adaptive paths

Detailed MC now records a six-bin histogram of the first triggering check and
the no-trigger shot count. Adding these counters does not consume random
numbers or change failure sampling. For each data qubit, report:

- Sum of idle durations in the circuits actually executed, in seconds.
- Sum of the individual idle-location Z probabilities. This is the expected
  number of physical idle-Z faults on that qubit, not a net Z probability,
  probability of any fault, or a logical error probability.

`idle_exposure` lists all seven possible paths: six triggered flagged prefixes
followed by the bare round, plus the full flagged round without a trigger.
The main per-round exposure is the path-frequency-weighted MC mean, with a
standard error from the sampled path distribution. The no-trigger and maximum
executed-path exposures are retained, rather than counting all potential
locations as if every circuit always ran. The report also gives readout-only
waiting seconds for each path (identical across data ions, since none is
measured), connecting the timing sweep to data idle exposure. The previous
exact boundary-subdivision convention remains unchanged.

### Artifacts and caching

`results/uniform_mc_calibration.json` and
`results/measurement_duration_sweep.json` checkpoint the expensive experiments.
Caches are accepted only when simulation-source fingerprint, parameters,
shot count, seeds and study grid match; use a new `--cache-dir` for a new
configuration. The final `results/ion_noise_results.json` contains all inputs,
checks, MC estimates, attribution, p_eff and ratios with confidence intervals,
pseudothreshold, exposure and the approved Stage 2 T2 sweep as a marked
reference. Original Stage 0–2 result files remain unchanged.

`figures/ion_noise_breakdown.png` groups raw P2 credit by location type: split a
mixed-type pair equally between its two types, and give a same-type pair full
credit to that type. These credits sum to P2; this is a stated attribution
convention. The figure labels P2 as secondary. The separate
`figures/ion_noise_calibration.png` shows the MC grid, p_L=p crossing, primary
targets and their sampling intervals. Stage 4 transport hooks are not added.
