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
targets and their sampling intervals. Stage 4 transport hooks are described below.

## Stage 4: optional transport hooks, synthetic validation only

`IdleIonNoise(params, transport_events=None, heating_penalty=None)` accepts the
scheduler interface. Defaults remain **no events** and **zero heating penalty**.
No trap geometry, device mapping, transport policy or hardware schedule is
created. The existing Stage 3 analysis remains the default CLI action.
Stage 4 is opt-in with `--stage 4`. H2 remains Stage 1 only because it lacks
sourced timing. Transport is available programmatically for a preset with all
required sourced timing inputs.

### Event convention and waiting time

Each event is `(stage, circuit_index, before_op_index, qubit, kind)`:

- stage is F or B; circuit_index is the original ROUND index, 0 through 5.
- before_op_index indexes the **original** `stabilizer_circuit` list before
  adding IDLE markers. This keeps scheduler inputs stable when timing changes.
  Index `len(original_ops)` means after that circuit's last op, before any
  successor circuit/adaptive transition.
- qubit is the moved ion, 0 through 8. kind is shuttle, split, merge or junction.
  Durations come exclusively from the corresponding preset fields.

Events at one insertion point execute sequentially in input order. The order
between distinct insertion points follows actual circuit execution. Every live
ion waits for each event's whole duration, **including the moved ion if live**.
Data are always live; ancillas are live after reset through their terminal
measurement. Events before an ancilla reset or after its measurement therefore
charge the live data/other ions, not that inactive ancilla. Transport itself has
no ideal logical action or separate X/Y/leakage channel.

The event duration is added to each live ion's pending idle gap and lumped into
the same Z channel before its next operation (or the circuit-end data tail).
Existing boundary subdivision remains exact. Event time is not added as an
extra noisy location per ion per transport action; it extends the appropriate
existing gap or creates one where the old gap was zero. Physical gate order
and the adaptive stopping rule are unchanged. Skipped circuits contribute
neither event time nor heating to executed gates.

### Next-executed-CX heating hook

`heating_penalty[kind]` is an additive **probability increment**, not a measured
energy or a calibrated fidelity model. Defaults and synthetic study values live
in `ion_params.py` with explicit placeholder comments. At an event, add that
increment to the moved ion's pending heating. At its next **executed** CX, add
all pending increments from both participant ions to the nominal p_2q and clear
both pending accumulators. Each event is consumed by only that CX. A reset does
not clear motional heating; no reset-induced cooling is assumed. Increments
that find no later CX in this one-round simulation have no realized gate
penalty, while their event duration is still charged. No cross-round motional
state is modeled.

Multiple events sum; moving both participants sums their increments on their
shared next CX. Negative/nonfinite increments, unknown event kinds, invalid
indices/ions and missing durations are rejected. If any reachable CX rate
would exceed one, construction fails rather than clipping it.

Because stopping early can send pending heating into B0 instead of the next
flagged check, B CX channels can depend on which flagged check triggered.
The model builds seven rate tables: each possible flagged stop and a no-trigger
reference. Flagged rates depend only on the preceding flagged prefix. The
no-trigger reference includes a virtual bare table after all flagged checks for
latent-location bookkeeping; those B gates never execute in no-trigger shots.

For branch-dependent cases, the MC dispatcher samples flagged faults, finds
the actual first trigger without applying a correction, then samples bare
faults using that branch's rate table and executes the existing EC simulator.
The preliminary frame pass consumes no RNG and changes no circuit/protocol.
A signal-provenance trigger split still uses reset/measurement/2q/idle/mixed:
transport waiting is idle, and the increased CX channel is 2q.

Pair counting chooses B payload probabilities from the branch produced by that
fault pair. It also uses that branch's nonidentity rates in the exact survival
product over every other potential location. The mean-rate approximation is
weighted per pair/branch. Bare-only singles remain skipped on noiseless input.
The support for fault testing is the union across possible branch channels.
The conditional pair counter and sampler are checked against an exact sparse
oracle with a forced flagged trigger and a branch-specific heated bare gate.

If all B rate tables coincide, the model uses the original independent-channel
counting/sampling kernels; the conditional sampler is also tested to give
identical seeded metrics in this limit. The symmetric synthetic fixture below
has this property, although asymmetric events exercise the full conditional
implementation in tests.

### Empty input regression and approved caches

No-event construction delegates to the original timing/channel path. The
independent MC/counting bodies are retained unchanged behind opt-in dispatchers.
With an empty list, circuit maps, ordered payload probabilities, P1/P2 and
breakdowns, seeded MC failure/trigger metrics and all existing results reproduce
Stage 3 exactly, even if a nonzero unused penalty map is supplied.

Stage 3 caches remain usable through an explicit frozen semantic signature
(`ion_cache.py`). It checks the original kernels, default dispatchers, protocol
and timing primitives, ordered baseline channels/circuits and study inputs
against the approved Stage 3 manifest before accepting the legacy signature.
Other changes fall back to the ordinary source fingerprint. No old cached
numbers or fingerprints are rewritten. Rebuilding the Stage 3 report with no
transport is checked byte-for-byte, including p_eff confidence intervals.

### Synthetic fixture and outputs

`--stage 4 --synthetic-transport` explicitly selects a test list of two syndrome-
ancilla shuttle events per flagged/bare check: before its first CX and before
its syndrome readout. There are 24 potential events; only the executed prefix
and, when triggered, bare events incur time. These are **synthetic markers**,
not a feasible routing plan or a hardware schedule. All use the sourced Jones
shuttle duration; the three penalty levels 0, 0.0005 and 0.001 are artificial
sensitivity values. They check strict monotonicity of raw P2, not a heating-law
fit or a monotonicity claim about a finite-shot estimator.

The checkpoint runner first replays empty transport at the approved nominal
Stage 3 seed, then samples each synthetic case with its own seed. MC remains
primary. Every case reports P1, raw/survival-corrected P2 and type/stage
breakdowns, MC p_L and stderr, trigger splits, failure share through bare,
separate flagged/bare-only/incoming fault counts, heated CX rates and adaptive
idle exposure including event time. Results are written separately to
`results/stage4_transport_results.json`. H2 and the nominal p_eff/ratios are
carried unchanged from Stage 3. The synthetic hook test is not assigned a new
hardware p_eff, and no extrapolation of the Stage 3 calibration is made.

Custom events may be supplied with `--transport-events events.json`, plus an
optional `--heating-penalty` JSON map. Stage 4 without either event input has
an empty event list; it does not automatically select the fixture. The original
Stage 0–3 numerical artifacts and figures remain unchanged.


## Stage A1: optional flag order and reset/readout layers

`FlagDesign` is an opt-in experimental view around the approved instrument.
`round_order` must be a permutation of `ROUND`, including the original support
indices. Both rounds execute that order; bare syndrome bits are stored by
support index, not appended by position. The triggering circuit's actual kind
and support select the original flagged table. No circuit, decoder, correction
rule or channel convention changes. Default order and SEQ delegate to the
original functions and retain the seeded random stream exactly. Other orders
reuse the independent or conditional counting/sampling kernel with private
function bindings; no module globals are mutated.

`SEQ` retains the original sequential timing. **PAR_MR is a modeling
assumption**, not a hardware schedule: greedily group adjacent reset-family
operations (R/RX), or adjacent measurement-family operations (M/MX), on distinct
ions. Each operation starts at the layer start, retains its original composite
duration, and the layer lasts the maximum duration. Gates remain sequential.
An RX still includes the reset and rotation channels; an MX still includes its
rotation and measurement channel. Only live-qubit waiting and resulting idle Z
locations change. A shorter reset's ion waits until its next gate; ancillas
cease to be live at their own measurement end. Data are live throughout.
Readout exposure uses the union of actual readout windows within a layer,
excluding an MX rotation from its own readout window.

For nominal Jones–Murali, each flagged check takes 1150 us under SEQ and
700 us under PAR_MR: max(50,55) + 6*40 + max(400,405) us. Bare Z/X checks
remain 610/620 us because each has only one reset and one measurement.
PAR_MR therefore differs from globally reducing the measurement-duration
parameter in the Stage 3 sensitivity sweep. Boundary subdivision and data
exposure are recomputed for each circuit and actual executed prefix. Intervals
below 1e-15 seconds are discarded solely to remove floating-point roundoff.

Transport events use circuit **position in the executed order**, with original
operation indices before idle insertion. An event at an insertion point is a
serial barrier between layers; terminal events are allowed. The same synthetic
waiting and prefix-dependent next-CX heating rules apply. Custom circuit maps
and idle inserters extend `TransportPlan` without changing its default path.
No hardware transport schedule is inferred.

The design grids are centralized in `ion_params.py`: UniformNoise at
{0.0005,0.001,0.0015}, and TWO_QUBIT_ONLY at {0.0005,0.001,0.002} as requested
experiment inputs, not new hardware parameters. TWO_QUBIT_ONLY uses the
existing CX depolarizing channel, with reset, readout and idle errors zero.
H2_LIKE stays timing-free, Stage 1 only; PAR_MR is rejected for untimed inputs.
Jones–Murali uses nominal sourced parameters and empty transport for both
reference rows. No physical placeholders were added.

A1 validates four orders under all nine preset/policy cases (36 rows). Exact
fault-support geometries can be reused across probability values: the saved
signature includes the order, operations, durations and every positive-probability
payload. There are 20 distinct exhaustive checks. Bare-only injected faults are
reported separately; noiseless input alone never executes the bare round.
The cross-check also forces a flagged trigger to execute each bare fault.

A1 MC is limited to the required baseline-order SEQ/PAR_MR timing references,
5M shots each with independent recorded seeds. MC is primary; raw P2 and exact
survival-corrected P2 are secondary attribution. Survival correction uses the
existing product over all other potential locations, including unexecuted
locations, and is not claimed to replace the adaptive MC estimate. JSON retains
both the exact correction and mean-probability approximation, cause-separated
trigger rates, trigger-position counts, failure share through bare, and all
actual-path idle-exposure tables. New files use the `stageA1_` prefix; all 16
approved result/figure files are SHA-256 checked against the pre-A1 commit.
No order-effect conclusion or new p_eff calibration is made until A2.

## Stage A2: order comparison and MC measurement

Stage A2 evaluates `ZZZXXX`, `XXXZZZ`, `ZXZXZX`, and `XZXZXZ`, preserving support
order within each kind. Every order uses both flagged and bare rounds in that
order. The nine comparison families are the three requested UniformNoise rates,
three requested TWO_QUBIT_ONLY rates, H2_LIKE (Stage 1 only, timing-free), and
nominal Jones–Murali separately under SEQ and PAR_MR. Transport is empty. Each
order is compared only to ZZZXXX in its own family; the timing assumption is
never changed within an order comparison.

### Sparse batched sampling

`design_mc.py` supplies a new sampler for this experiment, without modifying the
approved scalar instrument. For n shots at a location with total fault
probability p, draw K~Binomial(n,p), select a uniformly random K-subset of shot
indices, and draw each payload from the original normalized channel weights.
For any binary fault pattern with k hits the resulting probability is
P(K=k)/choose(n,k) = p^k(1-p)^(n-k), exactly the independent Bernoulli law.
Locations draw independently. Only shots with at least one fault need frame
propagation; omitted zero-fault shots remain in the denominator as successes.
Every nonzero shot executes the full adaptive protocol, including early exit,
no-trigger residual errors, flagged decoding and any triggered bare faults.
There is no pair-order truncation, rare-event approximation or change of channel.

Four linear component frames track reset, measurement, 2q and idle contributions
to the actual flagged trigger, with the same cancellation and mixed-cause rule
as the scalar instrument. Corrections use the original lookup tables. Syndrome
bit weights are the original (4,2,1) support addresses, independent of execution
position. Final logical-X-only, logical-Z-only and simultaneous XZ failures are
mutually exclusive and sum to p_L; inclusive X/Z rates may overlap through XZ.

Sampling uses its own recorded NumPy PCG64 streams, not the frozen scalar
random.Random streams. Each million-shot checkpoint stores raw accumulated
counts and RNG state atomically. Cache acceptance requires the same source,
ordered channels, timing policy, order, seed, NumPy version and chunk size.
Validation replays every positive payload and forced-trigger bare fault against
the scalar implementation across all 20 A1 support geometries, plus mixed
multi-fault cases and both-basis tableau checks. Sampling tests check location
independence, uniform and unequal payload weights, and exact cache resumption.

### Shot allocation, attribution and calibration

Sample allocation is fixed initially from deterministic P2, with variance margin:
per family n = max(5M, ceil_to_1M(4 / (0.02^2 * min_order P2))). P2 controls
allocation only; MC supplies every primary estimate. If measured combined stderr
exceeds 2% of the smaller of the two compared p_L values, the entire family is
extended in 1M-shot blocks. Differences use independent-stream combined stderr
sqrt(stderr_order^2 + stderr_baseline^2). The baseline compared with itself has
zero difference and zero comparison error. Significance means |delta|/stderr
>=3 for that comparison; no correction or claim of simultaneous family-wise
significance is made. Null results are retained.

Raw and exact survival-corrected P2 remain attribution only. For identical
uniform or CX-only support geometries across scan probabilities, exhaustive
counts can be reweighted exactly by p^2. Exact survival is (1-p)^(M-2) where M
is the number of active equal-rate locations (96 for uniform, 60 for CX-only);
zero-rate spectators contribute factors of one. The original all-location mean
approximation is saved separately. For H2/Jones and new orders, pairs are
explicitly enumerated. The two approved A1 Jones baseline pair tables are read
unchanged. Every row reports both P2/MC ratios and discrepancy in MC stderr.

The approved Stage 3 uniform MC grid remains the calibration instrument for all
orders. Seven lower points in `ion_params.py`, each with 5M shots, extend its
coverage; these are analysis scan inputs, not hardware parameters. Effective
rates use the approved piecewise log-log inverse and 20,000-draw independent
binomial bootstrap of targets and calibration, with the same calibration draws
shared across threshold ratios. Targets and bootstrap intervals must be
bracketed; extrapolation raises an error. P2-based p_eff = sqrt(P2/A_uniform)
is secondary. MC p_eff confidence intervals, p_eff/p_2q, and p_eff relative to
the uniform MC pseudothreshold are reported for every row.

A1's source identity excludes only the appended A2 analysis-grid section of
`ion_params.py`, retaining exactly its original bytes for all A1 inputs. The A2
fingerprint includes the entire parameter file and the new sampler. All 16
Stage 0–4 artifacts and four approved A1 JSONs are SHA-256 checked and untouched.
No decoder, check circuit, physical timing value or hardware schedule is changed.

## Stage A3: carried-frame memory and marginal error metrics

The user extends the earlier one-round scope to R={1,2,4,8,16}, flag orders
ZZZXXX/XXXZZZ/ZXZXZX, Uniform(0.001), H2 with timing-free Stage 1 channels,
and Jones–Murali under SEQ/PAR_MR. Each trajectory starts with ideal encoded
data and zero error frame. The exact output frame becomes the next round's
input. There is no re-preparation, ideal cleanup, frame normalization or added
waiting between rounds. Actual noisy protocol corrections occur each round;
one ideal decoder is used only for final observation. The same Pauli trajectory
evaluates both bases: logical X parity fails |0_L>, logical Z parity fails
|+_L>. Logical changes can cancel across rounds; first-passage failure is not
an absorbing logical-failure criterion.

### Reusable harness and sampling

A `RoundAdapter` supplies `initial_state(shots)` and `step(state,rng)` returning
`RoundStep(state,branch_taken,aborted,diagnostics)`. The `memory_trial` driver
carries this opaque state and calls a scorer once at the endpoint. An abort
can be recorded as an absorbing failure by a future retry/preparation adapter;
that adapter must retain its abort status. No Stage B circuit or verification
implementation is included in A3.

`FlagMemoryRound` reuses the approved circuits, channels and decoder tables.
Its sparse shortcut skips only shots with both zero incoming frame and no new
faults; carried errors are never discarded. Ancillas reset only through actual
circuit reset operations. Noise is sampled at every potential location in every
round; the unchanged adaptive protocol ignores latent unexecuted faults. This
is full finite-rate MC, with no pair-order truncation. Existing live-ion idle
insertion covers every executed check; there is no extra round-boundary time.

Four current-round components (reset, measurement, 2q, idle) and a fifth incoming
`carry` frame are replayed through the actual flagged prefix. Exclusive mixed
causes, overlapping contributors and fractional credits are saved. Historical
noise types behind a nonlinear earlier correction are not arbitrarily assigned
to current locations. Trigger rates are per executed round; any-bare probability
is per trajectory. Failure shares through any bare round and the last bare
round are both reported. Data idle exposure is averaged over actual round paths;
independent-round exposure stderr is omitted because paths within a trajectory
can correlate.

Every design/R uses an independent NumPy stream. One million trajectories per
point gives 60M endpoints and 372M noisy rounds. Atomic 250k-trajectory caches
store counts and RNG state. Identity includes the memory kernel, physical
parameters/channels, order, R, seed, NumPy version and chunk size. All 68 approved
Stage 0–A2 artifacts remain SHA-256 identical. Older source identities exclude
only validated literal later analysis lengths, never physical/runtime overrides.
Fit code has a separate source hash; the user's metric clarification requires
no new noise sampling.

### Primary marginal fits and combined metric

As clarified by the user, fit epsilon_X and epsilon_Z **separately** to
q_b(R)=[1-(1-2 epsilon_b)^R]/2 for R={2,4,8,16}, with binomial likelihood.
R=1 is reported but excluded. The union p_L is descriptive at every R and is
**not fitted**. Epsilon_tot is defined as

    epsilon_tot = 1 - (1-epsilon_X)(1-epsilon_Z).

This user-defined combined metric is used for every primary design delta and
sigma; it is not substituted for an observed union endpoint probability.
Stage B3's primary design metric is now epsilon_tot from the marginal memory
fits, with one-round union p_L secondary. Epsilon_X/Z and their fit quality
remain required diagnostics.

The two marginal fits share the same Pauli trajectories. Their covariance must
therefore be propagated. Let H_b be the observed negative log-likelihood
curvature and w_bR = q'_b(R)/[H_b q_b(R)(1-q_b(R))]. Independent R samples give

    Var(epsilon_b) = sum_R w_bR^2 N_R p_bR(1-p_bR)
    Cov(epsilon_X,epsilon_Z) = sum_R w_XR w_ZR N_R (p_XZR-p_XR p_ZR).

Here p_XZR is the observed probability of failure in both bases. The marginal
stderr uses this delta-method variance of the finite-window MLE; the expected
Fisher stderr is retained as a secondary check. Then

    Var(epsilon_tot) = (1-epsilon_Z)^2 Var(epsilon_X)
                    + (1-epsilon_X)^2 Var(epsilon_Z)
                    + 2(1-epsilon_X)(1-epsilon_Z) Cov(epsilon_X,epsilon_Z).

This retains positive X/Z sampling correlation instead of assuming independent
marginal estimators. Normal 95% intervals, covariance and correlation are saved.
The reported sampling error concerns the fitted finite-window metric; poor
fit quality still limits a stationary-rate interpretation. No one-round p_eff
calibration is repurposed as a memory calibration.

### Fit quality and the R=2 boundary check

Each marginal fit reports Pearson chi-square (3 degrees of freedom), its
p-value, binomial deviance and standardized residuals. An independent R≥4 fit
predicts the held-out R=2 point. The report gives the deviation in R=2 MC
error bars and also the combined uncertainty including the trend prediction.
One error bar means the reported one-sigma MC stderr. Every case outside that
bar has its R≥4 fit reported; all other R≥4 fits are shown too, avoiding
selective comparison windows. R≥8 is an additional diagnostic, not a replacement
for the requested R≥2 primary metric. Pairwise epsilon_tot and epsilon_Z
comparisons are supplied for both primary and later windows.

Dropping R=1 alone is not asserted to eliminate all residual-syndrome boundary
effects. Fits with poor goodness of fit are reported as such. A Pearson-scaled
stderr is only a heuristic diagnostic, not a validated bias correction. Primary
deltas use independent-design combined propagated sampling stderr and the
specified 3-sigma rule, without a simultaneous/family-wise claim. No infinite-R
extrapolation is established by this experiment.

### Exact one-round alternating-order mechanism

`pair_mechanism.py` evaluates all 431,280 distinct-location payload pairs under
Uniform(0.001) for both ZZZXXX and ZXZXZX. Physical identity is matched by stage,
kind/support, original operation index and payload, remapping executed positions.
Every class retains location/operation types, both position maps, both trigger
paths, and new/removed/common malignant weights. Full positive and negative
classes plus scalar propagation traces are saved. Raw P2 reproduces A2;
exact survival correction is (1-0.001)^94 under the approved full potential-
location convention. These one-round attributions are compared only against
one-round MC; they are not multiplied by R or used as epsilon estimates.

## Stage B1: published verified preparation (standalone)

### Source, code convention and circuit

The preparation is transcribed from [Goto, Scientific Reports 6, 19578 (2016),
Fig. 1(c)](https://doi.org/10.1038/srep19578). The prompt's arXiv:1610.03309 is
an unrelated hybrid-computation paper. The publisher's
[2022 correction](https://doi.org/10.1038/s41598-022-06633-6) fixes the old PDF
footer's volume number; it does not change the circuit. The construction and
Hadamard-dual CSS principle are also described in
[Peham et al., PRX Quantum 6, 020330 (2025), Example 3/Fig. 4](https://doi.org/10.1103/PRXQuantum.6.020330).

Use local block labels 0..6, physically mapped to 7..13, with verifier 14.
The eight-CNOT encoder is

    (1,0), (3,5), (2,6), (1,4), (2,0), (3,6), (1,5), (6,4).

Initial local inputs are |0,+,+,+,0,0,0>. The verifier is reset to |0>.
CNOTs (0,V),(5,V),(6,V), then M(V), check logical Z0 Z5 Z6. Accept only the
reported bit 0; no correction is applied. Disjoint gates drawn in one paper
column are linearized top-to-bottom and remain serial in the timing model.
All eight resets occur first. This program-order convention was fixed before
single-fault validation and is not a hardware schedule.

Propagating initial X generators on local 1,2,3 produces supports {0,1,4,5},
{0,2,4,6}, {3,4,5,6}; their span equals the frozen X-stabilizer span. Z0 Z5 Z6
commutes with all code generators and has odd weight, so it is a logical Z
representative. The generated superposition has eigenvalue +1 for it.
The unverified encoder's dangerous X components fall into two bad cosets,
represented by X4 X6 (equivalently X0 X2) and X1 X5 (equivalently X0 X4).
Both anticommute with this verification observable. Independent exhaustive
propagation finds 30 encoder-payload instances in these cosets. The verifier
therefore rejects their one-fault instances. Its own CNOT faults can add at
most one X on the output block.

The |+_L> preparation is the full Hadamard conjugate: interchange R/RX and
M/MX, reverse every CNOT. It checks X0 X5 X6 with V in |+> and MX readout.
This uses existing reset/rotation primitives, not an extra noiseless or
unsourced H gate. Dangerous Z replaces dangerous X. Physical channel weights
retain their original basis dependence; acceptance is not assumed equal for
ion presets or Z-only idle noise.

### What the single-fault guarantee means

Reduce the dangerous component modulo the **original code stabilizers only**;
do not quotient out the dangerous logical operator. Accepted single-fault
outputs have dangerous reduced weight <=1. Frame/tableau checks compare the
six code-stabilizer signs, prepared logical-observable sign, verification
outcome and a dangerous-weight classification independently reconstructed from
the tableau signature. No correction table is trained on failing faults.

The other component reaches raw and code-reduced weight 3 in the canonical
propagated Pauli frame. However, the prepared-state stabilizer additionally
contains logical Z for |0_L>, or logical X for |+_L>. Modulo that group the
orthogonal component has weight at most one. This is state degeneracy, not a
change of code or decoder. Z on an ancilla used as CX control does not spread
to the data target; X on an ancilla used as target does not spread to the data
control. B1 proves the stated preparation component criterion; the full noisy
EC round remains a separate B2 correctness test.

### Timing and scope

There are 8 resets, 8 encoder CX, 3 verification CX and 1 verifier readout.
SEQ retains serial operations. PAR_MR puts the eight adjacent disjoint resets
in one layer; gates remain serial, and there is only one measurement in this
gadget. Reset/readout operation durations and Pauli channels are exactly the
existing conventions, including RX/MX rotation contributions.

Only prepared/live ions receive idle Z noise. Block ions remain live through
the full verifier readout; the verifier ceases being live at its terminal
measurement. Idle waits split at reset/encoder/verification/readout boundaries
for attribution, without added time. Nominal Jones durations per attempt are
1.255 ms / 1.270 ms for zero/plus under SEQ and 0.895 ms / 0.900 ms under PAR_MR.
H2 stays timing-free. Transport is empty.

B1 is standalone: data 0..6 are reserved but untouched, and no data-storage
noise is simulated. A later B2 embedding must add data waiting throughout
every attempt, including rejected attempts. Reusing the standalone acceptance
or duration as a full EC memory result would omit that cost. All rejected
attempts are real operations and time; their count is included in standalone
retry resource statistics.

### Acceptance, Pauli simulation and attribution

The verifier's ideal outcome is deterministic. All decisions in B1 depend
only on that reported bit, never the ideal residual weight. Diagnostic block
readouts are random codewords: scalar/tableau tests compare H*m, not the raw
bitstrings. No ancilla/data coupling or EC correction is implemented in B1.

Pauli faults propagate linearly through this fixed Clifford attempt. Compile
each location/payload into its output X/Z frame and reported-verification-bit
influence. Finite-rate MC samples the unchanged independent channels and XORs
those influences. Full output frames can be returned for future reuse.
A 256-state XOR convolution over the dangerous 7-bit mask plus verification bit
provides an exact finite-rate benchmark, not an MC replacement. The parity
product for rejection gives a second acceptance check.

P1/P2 refer to the **unconditional per-attempt probability of an accepted
preparation with dangerous reduced weight >1**, not full-EC p_L. Exact survival
correction multiplies each pair by the product of (1-p_loc) over every other
prep location; the mean-rate approximation is also retained. Both are compared
to the corresponding MC preparation-hazard rate. For example, extra orthogonal
idle faults may be irrelevant to this hazard but are excluded from an
exactly-two-fault probability; a survival-corrected P2 is not the full
finite-rate hazard. No memory epsilon or p_eff is inferred from these numbers.

Verification rejection is the B1 branch-trigger metric, split by reset,
measurement, 2q, idle and mixed causes using linear signal contributions.
Flag-trigger and bare-failure-share metrics are not applicable to standalone
preparation. The original low-count H2 hazard samples are retained; fixed,
independent 20M-attempt checks are additionally reported, rather than silently
replacing those samples.

`sample_retries` requires an explicit **total-attempt** cap. B1 reports both
three-total-attempt and initial-plus-three-retry (four-total-attempt) diagnostics;
the ambiguous “max_retries=3” convention is left explicit for review before
B2. A failed verification discards and resets the whole block; a dangerous
multi-fault state with verifier 0 is still accepted. Exhaustion is counted as
a standalone failure. For per-attempt rejection q, exhaustion is q^K, mean
attempts is sum(q^i, i=0..K-1), and the selected dangerous-state probability is
the one-attempt accepted hazard times that same sum. Finite-rate and MC values
are both saved. These standalone diagnostics do not simulate idle evolution
of an external data block.

### Later B3 decisions (recorded only)

The user selects XXXZZZ as the “flag, best order” arm under **every** preset,
including nonsignificant cases. Both Steane half-orders must be run. In the
prompt's definitions, the X-error-syndrome half uses |+_L> as target with Z
readout; the Z-error-syndrome half uses |0_L> as control with X readout. Primary
comparisons use A3 epsilon_tot, always with epsilon_X and epsilon_Z. No B3 arm
comparison, half-order round assembly, or later-stage sampling occurs in B1.
