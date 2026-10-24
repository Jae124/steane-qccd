# steane-qccd

Steane [[7,1,3]] fault-tolerant gadgets, their logical error rates, and (next) their
schedules on a QCCD trapped-ion architecture. Master's thesis work, Kubiatowicz group,
UC Berkeley.

## Layout

```
src/
  steane_ft.py              Layer 1: flagged EC round (Chao–Reichardt) + single-fault test
  steane_noise.py           Layer 2: per-location noise model, fault-pair counting, Monte Carlo
  ion_params.py             Centralized trapped-ion presets and fidelity conversion
  ion_noise.py              Stages 1–4: channels, checks, Monte Carlo, p_eff and transport input
  ion_timing.py             Sequential timing, idle Z channels and adaptive exposure
  ion_calibration.py        Log-log inversion and joint binomial uncertainty
  ion_analysis.py           Stage 3 calibration, readout sweep, results and figures
  ion_transport.py          Optional event hooks and synthetic validation
  ion_cache.py              Approved baseline cache compatibility checks
  plot_ec_threshold.py      p_L vs p figure
  draw_flagged_circuit.py   Clean drawing of the flagged Z4Z5Z6Z7 circuit
results/                    Generated data (.json, .stim)
figures/                    Generated figures (.png, .svg)
```

## Run

Legacy uniform/figure scripts run from `src/` and write into `results/` and
`figures/`. The ion CLI below also runs from the repository root.

```
pip install -r requirements.txt
cd src
python steane_ft.py            # expect: 0 / 1218 failures with flags, 48 / 1218 without
python steane_noise.py         # expect: A ≈ 768, Monte Carlo agrees at low p (~1 min)
python plot_ec_threshold.py
python draw_flagged_circuit.py
```

## Current results (one EC round, uniform circuit-level noise)

- Fault-tolerant: no single fault causes a logical error.
- p_L ≈ 768 p² (exact leading order over 96 fault locations); pseudothreshold ≈ 1.3 × 10⁻³.
- 90% of A comes from fault pairs that hit the unflagged follow-up round.
- By operation type: 389 p₂q² + 178 p₂q·p_meas + 178 p₂q·p_reset + 24 (meas/reset pairs).

## Not modeled yet

Shuttling / split / merge error and noisy logical encoding. Uniform
noise remains the default; operation-specific ion noise is opt-in. Pseudothreshold not yet checked against Chao–Reichardt.

## Next

1. Sourced transport/heating inputs and scheduler integration (no hardware schedule yet).
2. Other gadgets: Goto-verified |0_L⟩, transversal H / S / CNOT, logical measurement.
3. QCCD schedule per gadget → idle and shuttle error per location.

## General channels (Stage 0 checkpoint)

Noise models can now expose `channel(loc)`, a dictionary of fault payloads and
probabilities. `UniformNoise(p)` reproduces the original uniform model;
`NoiseModel(p)` remains compatible. `count_faults(noise)` reports P1, P2, and
logical/type/stage breakdowns; `monte_carlo(noise, shots, rng)` samples those
same channels. See [docs/noise_model.md](docs/noise_model.md) for definitions and
[docs/STAGE_LOG.md](docs/STAGE_LOG.md) for checkpoint results.

From the repository root:

```sh
.venv/bin/python -m pytest -m "not slow"  # fast correctness tests
.venv/bin/python -m pytest               # includes exhaustive pairs + 3M shots
```

The full suite includes the exact seeded baseline: 542 failures in 3,000,000
shots at p = 0.0005, p_L = 0.000180667 ± 0.000007760. The second-order estimate
is P2 = 0.000191978, with P1 = 0 and A = 767.911111. The slow checks took about 64 seconds on the checkpoint machine
(and may take longer on other machines).


## Trapped-ion operation rates (Stage 1 checkpoint)

```sh
# From the repository root; defaults to both presets, 3M shots each, seed 2026.
.venv/bin/python src/ion_noise.py --stage 1
# Or select one preset (also works as `python ion_noise.py` from src/):
.venv/bin/python src/ion_noise.py --stage 1 --preset H2_LIKE --output results/h2_stage1.json
.venv/bin/python src/ion_noise.py --stage 1 --preset JONES_MURALI_2026 --output results/jones_stage1.json
```

The CLI verifies channel fault tolerance and frame/tableau agreement before
counting pairs and sampling. With `--stage 1` it writes `results/stage1_ion_noise.json` by default.
At Stage 1 it models operation errors only; idle/transport and p_eff belong to later
checkpoints. Preset sources and conventions are in
[docs/noise_model.md](docs/noise_model.md). Jones–Murali gate probabilities are
placeholders; the H2-like depolarizing convention remains unverified.

| Preset | P1 | Raw P2 | Survival-corrected P2 | Monte Carlo p_L ± standard error |
| --- | ---: | ---: | ---: | ---: |
| JONES_MURALI_2026 | 0 | 0.00167173 | 0.001419706 | 0.001592 ± 2.3e-05 |
| H2_LIKE | 0 | 0.000884656 | 0.000805656 | 0.00087 ± 1.7e-05 |

The survival correction uses each pair's exact product over every other
potential location. It counts exactly two effective fault assignments across
the full pool of 96 locations; Monte Carlo includes additional assignments,
including assignments in skipped circuits. Consequently the corrected value is
a subset contribution, not the full finite-rate prediction. See the checkpoint
log for both comparisons and the operation-type breakdown.

## Idle dephasing (Stage 2 checkpoint)

```sh
# Stage 2: Jones–Murali Stage 1 reference + idle model + T2 sweep,
# and H2-like Stage 1 only. Monte Carlo is the primary p_L estimate.
.venv/bin/python src/ion_noise.py --stage 2
.venv/bin/python src/ion_noise.py --stage 2 --preset JONES_MURALI_2026
# Explicitly reproduce timing-free Stage 1 with the new trigger metrics:
.venv/bin/python src/ion_noise.py --stage 1 --preset H2_LIKE --output results/h2_stage1_metrics.json
```

`src/ion_timing.py` inserts Z-only idle channels for sequential operation timing.
Gaps are lumped within each circuit and subdivided at adaptive boundaries; both
sides of a cross-boundary wait are charged only when their circuit executes.
This preserves the exact dephasing channel. See
[the timing assumptions](docs/noise_model.md#stage-2-sequential-idle-dephasing-and-reporting-decisions)
for the boundary convention and its effect on P2 attribution.

Results go to `results/stage2_ion_noise.json`. Nominal Monte Carlo runs use
3,000,000 shots each; other T2 points use 1,000,000 shots each (override with
`--shots` and `--sweep-shots`). The file includes raw and survival-corrected P2,
the trigger rate and category splits, and the share of failures that execute
the bare round. H2-like stays timing-free and is explicitly marked “Stage 1
only.” The updated Stage 0 uniform trigger metrics are stored separately in
`results/stage0_trigger_metrics.json`.

The Stage 3 analysis below defines p_eff primarily from Monte Carlo; P2-based
p_eff is a secondary check.

Nominal Jones–Murali Stage 2 result (T2 = 2.2 s):

- Monte Carlo p_L = 0.002007 ± 2.58e-05 (3M shots).
- Trigger rate = 10.50293%; 98.7211% of failures went through the bare round.
- Attribution: raw P2 = 0.0021817915, survival-corrected P2 = 0.0018200904.
- Flagged single faults: 0/1332; bare-only (skipped on noiseless input): 0/888; incoming errors: 0/42.
- Full metrics, trigger category splits, and T2 sweep: [Stage 2 report](docs/STAGE_LOG.md#stage-2--sequential-idle-z-noise-and-primary-monte-carlo-reporting).

## Effective uniform rate (Stage 3 checkpoint)

```sh
# Default stage is now 3. Runs/resumes the grid and readout study, then reports.
.venv/bin/python src/ion_noise.py
# Also usable from src/: python ion_noise.py --preset H2_LIKE
# New sampling configurations require a separate cache directory:
.venv/bin/python src/ion_noise.py --grid-shots 2000000 --cache-dir results/stage3_2M_grid
```

The uniform grid has 11 points at 1M shots each, with six-point neighborhoods
bracketing both preset targets. p_eff is obtained by log-log interpolation.
The uniform pseudothreshold solves uniform MC p_L=p on that curve. The reported
95% intervals propagate binomial sampling uncertainty from the target and the
uniform grid, including their shared calibration in threshold ratios. These
are sampling intervals conditional on the interpolation model.

Jones–Murali measurement-duration sensitivity uses 0/100/200/400 μs at 3M
shots each, keeping measurement error and other parameters nominal. It also
reports the nominal per-data-qubit idle seconds and summed local Z probability,
weighted by the actual adaptive-path histogram. H2 remains **Stage 1 only**;
its approved 3M-shot timing-free target is reused.

Outputs:

- `results/ion_noise_results.json`: both presets, p_eff and confidence intervals,
  p_eff/p_2q and p_eff/pseudothreshold, pair attribution, timing sensitivity,
  idle exposure, trigger metrics and the Stage 2 T2 sweep as reference.
- `results/uniform_mc_calibration.json`, `results/measurement_duration_sweep.json`:
  resumable raw experiments with seeds, parameters, checks and source fingerprints.
- `figures/ion_noise_breakdown.png`: secondary P2 credit by location type
  (mixed-type pairs split equally).
- `figures/ion_noise_calibration.png`: uniform MC calibration and primary targets.

For the two expensive parts independently and a report-only rebuild:

```sh
.venv/bin/python src/ion_analysis.py --part calibration
.venv/bin/python src/ion_analysis.py --part measurements
.venv/bin/python src/ion_analysis.py --part report
```

See [docs/noise_model.md](docs/noise_model.md#stage-3-effective-uniform-rate-from-monte-carlo)
for uncertainty, attribution and counterfactual-timing definitions. Stage 4
transport hooks are described below; no hardware schedule is supplied.

Stage 3 headline results (95% sampling intervals):

| Preset | MC p_eff [95% CI] | p_eff / p_2q | p_eff / uniform pseudothreshold |
| --- | ---: | ---: | ---: |
| JONES_MURALI_2026 | 0.001675669 [0.001638817566, 0.001710165456] | 1.67567 | 1.25529 |
| H2_LIKE (Stage 1 only) | 0.0010734088 [0.001032147376, 0.001135605607] | 0.838601 | 0.804118 |

Uniform pseudothreshold = 0.001334889, 95% CI [0.00129130961, 0.00145135925]. Reducing readout duration from 400 to 0 μs lowers Jones–Murali p_L by 17.82%, accounting for about 81.1% of the nominal idle-induced logical increment under this timing model. This changes waiting for all live ions, including ancillas. Full readout p_L/p_eff tables, trigger metrics, per-data exposure and assumptions are in [the Stage 3 report](docs/STAGE_LOG.md#stage-3--monte-carlo-effective-uniform-rates-and-readout-sensitivity).

## Transport hooks (Stage 4 checkpoint)

The optional input is a list of `(stage, circuit_index, before_op_index, qubit,
kind)` events, with kind `shuttle`, `split`, `merge` or `junction`. Indices refer
to the **original op lists before IDLE insertion**. Transport duration adds idle
Z noise on every live ion; optional heating increments affect the next executed
CX on the moved ion, including across adaptive/circuit boundaries. Defaults are
an empty event list and zero heating. H2 remains **Stage 1 only**.

```sh
# Explicitly run the synthetic test fixture and placeholder heating sweep:
.venv/bin/python src/ion_noise.py --stage 4 --synthetic-transport
# Empty transport regression; no automatic schedule is generated:
.venv/bin/python src/ion_noise.py --stage 4
# Future scheduler/user-supplied event list, not a schedule synthesized here:
.venv/bin/python src/ion_noise.py --stage 4 --transport-events events.json --heating-penalty '{"shuttle":0.001}'
```

Programmatic use:

```python
noise = IdleIonNoise(params, transport_events=events,
                     heating_penalty={"shuttle": 0.001})
counts = count_faults(noise)
metrics = monte_carlo(noise, shots, rng, return_details=True)
```

`src/ion_transport.py` validates inputs, builds waiting/heating tables and handles
branch-dependent channels. `src/ion_cache.py` preserves approved Stage 3 caches
only when the baseline semantic signature matches. The independent MC/counting
kernels, decoder and gate order remain unchanged.

The checkpoint report is `results/stage4_transport_results.json`. The fixture's
24 potential events and nonzero heating penalties are **synthetic interface
tests only**, not a hardware schedule or calibrated heating model. Details and
scope are in [docs/noise_model.md](docs/noise_model.md#stage-4-optional-transport-hooks-synthetic-validation-only).

Stage 4 synthetic checkpoint: no-event results exactly reproduce Stage 3,
including seeded MC and rebuilt p_eff confidence intervals. All 61 tests pass
(fast suite: 55 tests in 49.84 s). Synthetic flagged faults: 0/1344; bare-only
faults: 0/912 separately; incoming: 0/42; no simulator mismatches.

| Per-event shuttle penalty | MC p_L ± stderr (3M shots) | Raw P2 | Survival-corrected P2 |
| ---: | ---: | ---: | ---: |
| 0 | 0.00208233333 ± 2.63e-05 | 0.00218859849 | 0.00182533392 |
| 0.0005 | 0.00226233333 ± 2.74e-05 | 0.00246016453 | 0.00202844645 |
| 0.001 | 0.00256566667 ± 2.92e-05 | 0.00274780167 | 0.00224023113 |

Raw P2 increases strictly over the synthetic penalty sweep. Full trigger
splits, type/stage attribution, event list and actual-path exposure are in
[the Stage 4 report](docs/STAGE_LOG.md#stage-4--optional-transport-interface-and-synthetic-validation).


### Stage A1: optional flag check order and timing policy

`src/flag_design.py` provides `FlagDesign(noise, round_order=None,
timing_policy="SEQ")`. The order is a permutation of the six `(kind,
support_index)` entries in `steane_ft.ROUND`, used for both flagged and bare
rounds. Defaults delegate to the approved instrument. Individual check circuits,
lookup tables and the adaptive rule are unchanged.

```python
from flag_design import CHECK_ORDERS, FlagDesign, TWO_QUBIT_ONLY
from ion_noise import IdleIonNoise
from ion_params import JONES_MURALI_2026

design = FlagDesign(IdleIonNoise(JONES_MURALI_2026),
                    CHECK_ORDERS["ZXZXZX"], timing_policy="PAR_MR")
checks = design.check_faults()
```

`PAR_MR` is an explicit modeling assumption: adjacent resets or adjacent
measurements on distinct ions share a layer; gates remain sequential. It uses
only the existing sourced durations. H2 remains timing-free (Stage 1 only).
Transport indices identify positions in the **executed** check order and
original operation indices before idle insertion; transport remains synthetic.

Reproduce this checkpoint with:

```bash
.venv/bin/python src/check_order.py --part checks
.venv/bin/python src/check_order.py --part reference  # cached 5M-shot SEQ/PAR_MR arms
.venv/bin/python src/check_order.py --part report
.venv/bin/python -m pytest -q
```

The A1 report is [docs/STAGE_LOG.md](docs/STAGE_LOG.md#stage-a1--parameterized-flag-check-order-and-timing-policy),
with machine-readable results in `results/stageA1_validation.json`. All 36
order/preset/policy configurations pass the fault and simulator checks. A1
contains only baseline-order timing reference MC; the order experiment (A2)
and Steane-style extraction (B) await review.

### Stage A2: check-order experiment

```bash
.venv/bin/python src/order_experiment.py --part all
# Independent resumable parts: pairs, mc, calibration, report
.venv/bin/python -m pytest -q
```

The experiment covers all four orders under the three uniform rates, three
CX-only rates, timing-free H2, and **both SEQ and PAR_MR for Jones–Murali**:
36 rows. Every difference uses the baseline order under the same preset and
timing policy. The primary report is `results/stageA_check_order.json`; the
nine-panel figure is `figures/check_order_comparison.png`. Detailed conclusions
and checkpoint validation are appended to `docs/STAGE_LOG.md`.

`src/design_mc.py` batches independent Pauli-frame shots with NumPy. It samples
a binomial number of faults at each location, chooses their shot indices
uniformly without replacement, and draws payloads using the unchanged channel
weights. This is the same independent-location distribution as scalar MC;
its separate RNG stream does not replace frozen seeded regressions. Exhaustive
payload and mixed-fault replays are compared shot by shot with the scalar
simulator, and tableau tests cover both logical bases. Transport is empty.

`results/stageA2_cache/` stores pair counts and 1M-shot MC chunks, with RNG state,
seed, accumulated counts, channel/configuration hashes and simulation source
fingerprints. Interrupted runs resume from the last complete chunk. The approved
Stage 3 uniform calibration is read unchanged and extended below its original
range with seven 5M-shot points; p_eff uses log-log MC interpolation with joint
target/calibration bootstrap intervals and no extrapolation. Stage B is deferred.

### Stage A3: multi-round memory (primary metric for Stage B3)

Stage B3 must use **epsilon_tot = 1 - (1-epsilon_X)(1-epsilon_Z)** as its
primary design metric. Epsilon_X and epsilon_Z are fitted separately to the
two logical-basis marginals over R≥2. Propagated uncertainty includes their
shared-trajectory covariance. One-round union p_L is secondary; union endpoint
probabilities at every R are descriptive and are **not fitted**.

```bash
.venv/bin/python src/memory_experiment.py --part mc
.venv/bin/python src/pair_mechanism.py
MPLCONFIGDIR=/tmp/steane_matplotlib .venv/bin/python src/memory_experiment.py --part report
.venv/bin/python -m pytest -q
```

A3 chains R={1,2,4,8,16} noisy rounds, carrying the output Pauli frame without
intermediate ideal decoding or extra time. It covers ZZZXXX, XXXZZZ, ZXZXZX;
Uniform(0.001), timing-free H2, and Jones–Murali SEQ/PAR_MR. Results and exact
matched-pair mechanisms are in `results/stageA3_memory.json` and
`results/stageA3_pair_mechanism.json`; the checkpoint is in `docs/STAGE_LOG.md`.
`figures/memory_epsilon_tot.png` shows the primary metric. Three
`memory_rounds_*.png` figures show basis-fit curves and descriptive union data.

`src/memory_mc.py` provides the reusable `RoundAdapter`/`RoundStep` interface
and `memory_trial` driver. The adapter supplies one noisy round; the driver
carries opaque state and scores only at the final endpoint. Only the flag
adapter is implemented here. An absorbing abort hook is available for future
verified-preparation adapters. Current-round reset/readout/2q/idle trigger
signals are distinguished from incoming `carry`; any-bare and last-bare failure
shares are separate. Caches save 250k-trajectory chunks and RNG state atomically.

`src/memory_fit.py` reports marginal goodness of fit, an R=2 comparison with
the independently fitted R≥4 trend, and R≥4/R≥8 sensitivity fits. The R≥4 fits
are shown for all cases, including every R=2 point outside its one-sigma error
bar. Inspect goodness of fit before interpreting a fitted finite-window metric
as a stationary rate. Stage B has not started.

### Stage B1: standalone verified ancilla preparation

`src/verified_prep.py` implements Goto's eight-CNOT encoder followed by a
three-CNOT logical verification measurement, plus its Hadamard dual. It prepares
|0_L> with dangerous X errors checked, and |+_L> with dangerous Z errors checked.
The published circuit uses the same stabilizer row span as the frozen code;
no flag circuit or decoder is modified. The correct source is
[Goto, Scientific Reports 6, 19578 (2016), Fig. 1(c)](https://doi.org/10.1038/srep19578).
The arXiv identifier in the attached design prompt belongs to a different paper.

```bash
.venv/bin/python src/prep_experiment.py --part checks
.venv/bin/python src/prep_experiment.py --part mc --retry-caps 3 4
MPLCONFIGDIR=/tmp/steane_matplotlib .venv/bin/python src/prep_experiment.py --part report --retry-caps 3 4 --retry-cap-selection "Both total-attempt caps are standalone diagnostics; the max_retries convention remains for review before B2."
.venv/bin/python -m pytest -q
```

The report is `results/stageB1_preparation.json`, with the checkpoint in
`docs/STAGE_LOG.md` and an annotated circuit diagram in
`figures/verified_prep_circuits.png`. Each of 18 state/preset/policy rows has 5M
single-attempt MC samples. Exact finite-rate Pauli convolution is an independent
benchmark. Raw and survival-corrected P2 refer only to an accepted preparation's
dangerous-component error; they are not a full-EC p_L or memory epsilon.

The layout reserves data 0..6, ancilla block 7..13 and verifier 14. B1 executes
only the latter eight ions. Its timing includes block waiting during verifier
readout and every rejected retry. **Data waiting must be added when embedding
this gadget in B2.** `sample_attempts(..., return_arrays=True)` supplies full
output Pauli frames. Acceptance uses only the reported verification bit, never
an ideal error-weight test. Retry sampling takes an explicit total-attempt cap;
three- and four-attempt diagnostics expose the ambiguous “max_retries=3” wording.

Recorded for later B3, with no comparisons run in B1: include XXXZZZ as the
“flag, best order” arm under every preset, regardless of significance; run
Steane-style in both half-orders; use A3 epsilon_tot as primary and always
report epsilon_X and epsilon_Z. B2/B3 assembly and sampling have not started.
