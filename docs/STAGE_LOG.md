## Stage 0 — General per-location Pauli distributions

Changes:
- `src/steane_noise.py`: explicit validated channels, `UniformNoise`, compatible legacy `NoiseModel`, weighted P1/P2 and breakdowns, general Monte Carlo, cross-check return counts.
- `tests/test_noise_channels.py`: validation, simulator regressions, sparse unequal-weight channel, exhaustive pair and seeded Monte Carlo regressions.
- `pytest.ini`: source import path and slow marker.
- `requirements.txt`: pytest added.
- `docs/noise_model.md`: channel contract, weighted sums, sampling, assumptions.
- `results/stage0_channels.json`: exhaustive uniform-channel count at p = 5e-4.
- `README.md`: channel API, tests, checkpoint numbers.

Checks:
- Single-fault test: 0 / 1218 (flags on); 48 / 1218 flags ignored.
- Cross-check: 0 mismatches / 1176 runs.
- Uniform regression: exact legacy A dictionaries and 92840 malignant pairs across 96 locations.
- Single-fault enumeration in pair counter: 0 failures.
- Channel probabilities validated at every location; invalid distributions rejected.
- Sparse unequal-weight custom model checked against independently weighted explicit fault pairs and Monte Carlo.
- Fast suite: 7 passed in 7.75 s; slow suite: 2 passed in 64.11 s.
- P2 = A['any'] p^2 within 2e-12 relative tolerance (independent floating-point accumulation).

Numbers (p = 5e-4):
- Monte Carlo: 542 / 3,000,000 failures with `random.Random(2026)`; p_L = 0.000180666666666667 ± 0.000007759596772647 (binomial standard error), reproducing the stored result exactly.
- P1 = 0.
- P2 = 0.000191977777777423.
- A: any = 767.9111111103565, X = 385.06666666666246, Z = 447.0222222223297.
- P2 by type pair:
  - 2q*2q: 9.72111111109e-5
  - 2q*reset: 4.43833333333e-5
  - 2q*meas: 4.43833333333e-5
  - meas*reset: 3.0e-6
  - reset*reset: 1.5e-6
  - meas*meas: 1.5e-6
- P2 by stage pair: FB = 0.000173253333333044, FF = 0.000018724444444447.

Assumptions made:
- Location channels are independent; at most one payload occurs per location.
- P2 uses distinct-location pairs without no-fault survival factors, as specified.
- Keep legacy fault supports, adaptive execution, and ideal preparation/final decoding.
- Preserve the legacy RNG path for equal-weight channels; use conditional weighted sampling otherwise.
- Explicit zero-weight payloads remain enumerated; empty distributions skip enumeration.
- No physical parameters or hardware claims introduced in Stage 0.

Open questions for the human:
- None required to complete Stage 0. Stage 1 awaits checkpoint review.


## Stage 1 — Operation-specific rates and survival correction

Changes:
- `src/ion_params.py`: immutable named presets, documented provenance and placeholders, average-fidelity conversion.
- `src/ion_noise.py`: exact operation channels, exhaustive channel checks and diagnostics, Stage 1 CLI and Monte Carlo comparisons.
- `src/steane_ft.py`: optional multi-location injection mapping for both stages; existing single-fault API retained. No decoder, circuit order, or adaptive-rule changes.
- `src/steane_noise.py`: per-pair exact survival weighting, corrected breakdowns, and mean-rate approximation, available to all noise models.
- `tests/test_ion_noise.py`: channel composition, MX reduction against explicit Pauli injection, fidelity conversion, fault coverage and pair/survival checks.
- `tests/test_noise_channels.py`: uniform survival identity; remove test writes to tracked Stage 0 results.
- `results/stage1_ion_noise.json`: all parameters, checks, counts, breakdowns, seeded Monte Carlo and comparison statistics for both presets.
- `docs/noise_model.md`, `README.md`: channels, exact reductions, sources, reporting requirement and run instructions.

Checks:
- Per preset, single-fault test: 0 / 1908 channel fault runs (flags on), comprising 1152 flagged and 756 bare-only runs; incoming errors: 0 / 42 additional runs. Comparable flagged-plus-incoming count is 1194; the original exhaustive reset-Pauli superset still gives 0 / 1218.
- Per preset, cross-check: 0 mismatches / 1908 channel fault runs, plus 42 incoming trials and 756 forced-trigger bare-fault trials. Bare-only singles are skipped on noiseless input; the forced-trigger cases exercise actual bare execution.
- Cross-check compares all six residual stabilizer signs, the relevant logical sign, and ideal-decoder failure decisions.
- `fault_list` covers every nonzero payload emitted by both ion presets; no new op or payload needed.
- RX convolution includes cancellation; terminal MX reduction checked against explicit premeasurement Paulis with readout flips, both logical bases and incoming error.
- Uniform regression passes: 0 / 1218 flags on, 48 / 1218 flags ignored; 0 / 1176 frame/tableau mismatches; 96 locations, 92840 malignant pairs and all legacy A dictionaries exactly reproduced; seeded 3M-shot uniform Monte Carlo gives exactly 542 failures.
- Survival checks pass: sparse nonuniform pair-by-pair product, sum of corrected breakdowns, and uniform identity P2_survival = P2*(1-p)^94. Product implementation handles probability-one locations without division.
- Fast suite: 24 passed in 21.04 s. Slow suite: 4 passed in 79.29 s (includes the exact 3M-shot uniform Monte Carlo regression). Both 3M-shot ion Monte Carlo runs also completed via the CLI.
- `python steane_ft.py` run from `src/` reproduces the baseline encoder and 0/1218 versus 48/1218 checks. `compileall` and `git diff --check` pass.

Numbers:

| Preset | P1 | P2 | Exact survival-corrected P2 | Mean-rate approximation | MC p_L ± stderr (3M shots) |
| --- | ---: | ---: | ---: | ---: | ---: |
| JONES_MURALI_2026 | 0 | 0.00167172991184 | 0.00141970623299 | 0.00141774620854 | 0.001592 ± 2.3017859e-05 |
| H2_LIKE | 0 | 0.00088465601772 | 0.000805656034421 | 0.000805285053837 | 0.00087 ± 1.7021977e-05 |

Each preset uses an independent `random.Random(2026)` stream, 3,000,000 shots.

### JONES_MURALI_2026

Monte Carlo failures: 4776. Mean location probability: 0.001751555.

- Raw P2: estimate/MC = 1.050082; (estimate − MC)/stderr = +3.464.
- Exact survival correction: estimate/MC = 0.891775; (estimate − MC)/stderr = -7.485.
- Mean-rate correction: estimate/MC = 0.890544; (estimate − MC)/stderr = -7.570.

| Location type pair | P2 contribution | Survival-corrected contribution |
| --- | ---: | ---: |
| 2q*2q | 0.000388844444444 | 0.000329233033942 |
| 2q*meas | 0.000178139673778 | 0.000150830665731 |
| 2q*reset | 0.000888268146667 | 0.000755120444673 |
| meas*meas | 6.04005280053e-06 | 5.11412280417e-06 |
| meas*reset | 6.0239463472e-05 | 5.12099835986e-05 |
| reset*reset | 0.00015019813068 | 0.00012819798224 |

| Stage pair | P2 contribution | Survival-corrected contribution |
| --- | ---: | ---: |
| FB | 0.00157813896073 | 0.00134038358614 |
| FF | 9.35909511111e-05 | 7.93226468444e-05 |

### H2_LIKE

Monte Carlo failures: 2610. Mean location probability: 0.00099952748.

- Raw P2: estimate/MC = 1.016846; (estimate − MC)/stderr = +0.861.
- Exact survival correction: estimate/MC = 0.926041; (estimate − MC)/stderr = -3.780.
- Mean-rate correction: estimate/MC = 0.925615; (estimate − MC)/stderr = -3.802.

| Location type pair | P2 contribution | Survival-corrected contribution |
| --- | ---: | ---: |
| 2q*2q | 0.000637082737777 | 0.000580243667306 |
| 2q*meas | 0.000229493402396 | 0.000208962031121 |
| 2q*reset | 1.13447724692e-05 | 1.03199873373e-05 |
| meas*meas | 6.11688485249e-06 | 5.56813739076e-06 |
| meas*reset | 6.02859442279e-07 | 5.48255057233e-07 |
| reset*reset | 1.53607827272e-08 | 1.39562093085e-08 |

| Stage pair | P2 contribution | Survival-corrected contribution |
| --- | ---: | ---: |
| FB | 0.000770924081555 | 0.000702072815799 |
| FF | 0.000113731936165 | 0.000103583218622 |

Assumptions made:
- Independent local channels; one effective location each for RX and MX after exact local composition/marginalization. RX X outcomes remain explicit even though they stabilize |+>.
- MX rotation depolarization becomes an effective FLIP by anticommutation and XOR with readout error. The measured ion is discarded/reset before reuse. No extra ideal gates or changes to extraction are introduced.
- Effective CX has one two-qubit depolarizing channel; native MS/CX synthesis errors beyond the specified rates are not modeled.
- Jones–Murali gate rates are explicitly marked placeholders. H2-like gate rates are used directly, with the process-vs-average-infidelity convention unverified. Its readout rate is the source paper's compromise, not a verified hardware-wide scalar.
- H2 timing/T2 are unspecified in Stage 1; no timing is inferred from another preset.
- Exact survival factors include every potential location outside the selected pair, including skipped circuits. The corrected quantity counts exactly two assignments in that full latent pool; it is not the full Monte Carlo probability. Monte Carlo includes more assignments, some operationally irrelevant when their circuits are skipped.
- Raw P2 and the survival-corrected quantity are both compared with Monte Carlo without requiring either truncation to agree within sampling error. Here raw P2 is closer for both presets; both corrected values underestimate MC beyond sampling error.
- As requested in the Stage 1 authorization, every later checkpoint must retain both raw P2 and a survival-corrected value and compare both against Monte Carlo. This is implemented generically in the counter and documented in `docs/noise_model.md`.

Open questions for the human:
- H2 timing and an appropriate stochastic T2 will need an explicit assumption or source when Stage 2 starts; the cited H2-like paper supplies a coherent dephasing model rather than the requested T2 model.
- The H2 depolarizing convention and Jones–Murali heating-dependent gate fidelities remain unverified/placeholder as specified.
- Stage 2 has not started; stop here for review.


## Stage 2 — Sequential idle Z noise and primary Monte Carlo reporting

Scope: Stage 1 approved. Timing and idle noise enabled for Jones–Murali only;
H2-like is explicitly Stage 1 only. Stage 3 has not started.

Changes:
- `src/ion_timing.py`: sequential durations, live windows, lumped local waits, exact boundary subdivision, stable idle probability formula.
- `src/ion_params.py`: centralized requested T2 sensitivity sweep; H2 timing remains unspecified.
- `src/steane_ft.py`: IDLE execution/fault enumeration/Stim export, optional circuit map, trigger diagnostics; original circuit generation and decoder unchanged.
- `src/steane_noise.py`: timed circuit maps in frame/counting/Monte Carlo, idle type, optional detailed trigger and bare-failure metrics; legacy integer Monte Carlo API retained. Uniform CLI now reports the same metrics.
- `src/ion_noise.py`: Stage 2 CLI and T2 sweep, timing-free H2 marker, Stage 1 references with new metrics, separate flagged/bare-only/incoming fault totals.
- `tests/test_ion_timing.py`: time conservation, noiseless preservation, exact idle movement/composition, all idle faults, infinite-T2 regressions, trigger attribution including mixed/cancelled cases, reporting preserves seeded sampling.
- `results/stage2_ion_noise.json`: primary Monte Carlo estimates, pair attribution, sweep, trigger splits and bare-failure share for nominal/sweep/reference cases.
- `results/stage0_trigger_metrics.json`: Stage 0 uniform trigger backfill, separate fault counts and pair comparisons; previous approved result files left intact.
- `docs/noise_model.md`, `README.md`: timing and boundary assumptions, primary-estimate policy, trigger attribution convention and run instructions.

Checks:
- Jones–Murali flagged single faults: 0 / 1332, including 0 / 180 idle-Z runs (both logical bases).
- Bare-round-only single faults: 0 / 888, including 132 idle-Z runs; these locations do not execute on noiseless input and are excluded from the executed-single-fault headline.
- Incoming single data errors: 0 / 42.
- Frame/tableau cross-check: 0 mismatches / 2220 channel fault runs plus 42 incoming runs. Another 888 forced-trigger bare-fault trials exercise actual bare execution with 0 mismatches; these are two-fault agreement tests, not FT single-fault trials.
- All finite T2 sweep points have identical support and pass the same fault checks; P1 = 0 throughout.
- Noiseless IDLE-inserted run: no trigger, trivial syndromes, all six stabilizers and encoded logical sign preserved for both bases.
- Infinite T2: original location map returned, all Stage 1 pair/count values exactly match the approved JSON; seeded sampling (including new metrics) exactly equals Stage 1.
- Uniform regressions: original 0/1218 flags on, 48/1218 flags ignored, 0/1176 frame/tableau mismatches, exact A dictionaries/92840 malignant pairs, and exactly 542/3M failures at p=5e-4 pass.
- Full suite: 36 tests passed in 119.48 s; the subsequently added all-gap movement test passed separately. Final fast suite: 32 passed in 29.18 s (5 slow deselected). Total current tests covered: 37.
- `compileall` and `git diff --check` pass. Stage 1 Jones/H2 counts and seeded failure totals are independently checked against the approved results when generating this report.

Timing sanity:
- Measurement-only wait: 0.0004 s, p_idle = 9.09008269471e-05.
- Largest represented gap: 0.00115 s, p_idle = 0.000261295337314; location ['F', 0, 18, ['IDLE', 0, 0.00115]].
- A gap can span multiple operations. The largest gap is a full 1150 microsecond flagged circuit for an unused data ion, rather than a single 400 microsecond measurement. Its probability is correspondingly larger and consistent with the same T2 formula.
- 252 potential locations: 96 original operation locations plus 156 idle locations (90 flagged, 66 bare).

Numbers: Monte Carlo is primary; P2 values are attribution only. Standard errors are one binomial standard error.

| Preset/stage | Shots | Monte Carlo p_L ± stderr | Raw P2 | Exact survival-corrected P2 |
| --- | ---: | ---: | ---: | ---: |
| Uniform p=5e-4, Stage 0 | 3000000 | 0.000180666666667 ± 7.7595968e-06 | 0.000191977777777 | 0.000183161425518 |
| Jones–Murali, Stage 1 reference | 3000000 | 0.001592 ± 2.3017859e-05 | 0.00167172991184 | 0.00141970623299 |
| Jones–Murali, Stage 2 (T2=2.2 s) | 3000000 | 0.002007 ± 2.5839066e-05 | 0.00218179149881 | 0.00182009035186 |
| H2_LIKE, Stage 1 only | 3000000 | 0.00087 ± 1.7021977e-05 | 0.00088465601772 | 0.000805656034421 |

### Trigger metrics

All entries are probabilities per shot; the last column is conditional on logical failure.

| Preset/stage | Trigger rate ± stderr | Failures through bare / failures | Failure share through bare |
| --- | ---: | ---: | ---: |
| Uniform p=5e-4, Stage 0 | 0.0267566666667 ± 9.3167855e-05 | 541 / 542 | 0.99815498155 |
| Jones–Murali, Stage 1 reference | 0.0973543333333 ± 0.00017114951 | 4743 / 4776 | 0.993090452261 |
| Jones–Murali, Stage 2 (T2=2.2 s) | 0.105029333333 ± 0.00017701052 | 5944 / 6021 | 0.987211426673 |
| H2_LIKE, Stage 1 only | 0.0499503333333 ± 0.00012577135 | 2552 / 2610 | 0.977777777778 |

Exclusive trigger causes: categories are credited only when they contribute to an actually nonzero bit at the first trigger. Multicategory contributions go in mixed. Rows sum to the trigger rate.

| Preset/stage | Reset | Measurement | 2q | Idle | Mixed |
| --- | ---: | ---: | ---: | ---: | ---: |
| Uniform p=5e-4, Stage 0 | 0.005869 | 0.005966 | 0.014906 | 0 | 1.566666667e-05 |
| Jones–Murali, Stage 1 reference | 0.057184 | 0.01136966667 | 0.02860433333 | 0 | 0.0001963333333 |
| Jones–Murali, Stage 2 (T2=2.2 s) | 0.056923 | 0.011345 | 0.02852066667 | 0.007997333333 | 0.0002433333333 |
| H2_LIKE, Stage 1 only | 0.0005803333333 | 0.01172433333 | 0.03761433333 | 0 | 3.133333333e-05 |

The JSON also provides overlapping contributing-type rates and an equal-share allocation of mixed shots over the four types. This convention preserves four-way totals without claiming a unique cause for mixed events. RX/MX rotation contributions retain their reset/measurement grouping.

### Nominal Stage 2 pair attribution

P1 = 0. Raw and corrected pair sums are compared against the primary Monte Carlo estimate without requiring agreement at finite rates.

- Raw P2: 0.00218179149881; value/MC = 1.087091; (value − MC)/stderr = +6.765.
- Exact survival correction: 0.00182009035186; value/MC = 0.906871; (value − MC)/stderr = -7.234.
- Mean-rate correction: 0.00181512396895; value/MC = 0.904397; (value − MC)/stderr = -7.426.

| Location type pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| 2q*2q | 0.000388844444444 | 0.000323602954618 |
| 2q*idle | 0.000244229987957 | 0.00020308726939 |
| 2q*meas | 0.000178139673778 | 0.000148251372267 |
| 2q*reset | 0.000888268146667 | 0.0007422074391 |
| idle*idle | 6.2034119455e-05 | 5.15425354296e-05 |
| idle*meas | 3.40953488086e-05 | 2.83514271743e-05 |
| idle*reset | 0.000169702130754 | 0.000141680701929 |
| meas*meas | 6.04005280053e-06 | 5.02666828385e-06 |
| meas*reset | 6.0239463472e-05 | 5.03342626347e-05 |
| reset*reset | 0.00015019813068 | 0.00012600572103 |

| Stage pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| FB | 0.0019670291389 | 0.00164135880758 |
| FF | 0.000214762359919 | 0.000178731544274 |

### Jones–Murali T2 sensitivity sweep

Nominal point: 3M shots; other points: 1M shots each. Each uses `random.Random(2026)`. All P1 terms and single-fault failures are zero. Detailed attribution and trigger-cause splits for every point are in the JSON.

| T2 (s) | MC p_L ± stderr (primary) | Raw P2 | Survival-corrected P2 | Trigger rate | Failure share through bare |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 0.030785 ± 0.00017273472 | 0.0412733436375 | 0.0240682142065 | 0.252216 | 0.935293162254 |
| 0.5 | 0.004392 ± 6.6126472e-05 | 0.00484124242974 | 0.00380527399812 | 0.130454 | 0.97427140255 |
| 2.2 | 0.002007 ± 2.5839066e-05 | 0.00218179149881 | 0.00182009035186 | 0.105029333333 | 0.987211426673 |
| 10.0 | 0.001611 ± 4.0104921e-05 | 0.00177331324753 | 0.0015000552929 | 0.099156 | 0.98696461825 |

Assumptions made:
- Sequential operation timing, no transport. Data live throughout; ancilla waiting only between reset and measurement. RX/MX include their rotation duration as active time for that ion.
- Boundary subdivision explicitly cuts a maximal cross-circuit data gap at possible termination points. Within each circuit each untouched wait is lumped once. The dephasing semigroup identity and commutation make the resulting Monte Carlo channel exact; P2 attribution depends on this effective-location representation.
- All 252 potential locations enter each survival product outside the pair, including inactive locations; corrected P2 is a two-assignment subset contribution, not the primary p_L estimate. More faults (including faults on skipped circuits) explain why Monte Carlo need not equal it.
- H2_LIKE remains Stage 1 only with timing fields None; no H2 T2 sweep is invented.
- Trigger categories describe observed signal provenance. Mixed shots have no forced unique cause; equal sharing is reported as an optional four-category convention alongside the exclusive mixed category.
- Monte Carlo is primary from this stage onward. Stage 3 p_eff must be computed from Monte Carlo; P2-derived p_eff is only a secondary check. Stage 3 work has not begun.
- Published input rates and the existing gate-fidelity placeholders/convention caveats are unchanged.

Open questions for the human:
- No decision blocks completion of Stage 2. H2 timing stays deferred until a sourced preset is supplied.
- Boundary subdivision and mixed-trigger allocation are explicit modeling/reporting conventions for checkpoint review.
- Stop here for Stage 2 review.


## Stage 3 — Monte Carlo effective uniform rates and readout sensitivity

Scope: Stage 2 approved. Stage 3 only. Monte Carlo is primary, H2 remains
Stage 1 only for physical noise, and no Stage 4 transport work has started.

Changes:
- `src/ion_calibration.py`: log-log inverse calibration, monotone bootstrap draws, uniform MC pseudothreshold and propagated target/grid confidence intervals, including covariance in threshold ratios.
- `src/ion_analysis.py`: cached/resumable calibration and readout experiments, complete Stage 3 JSON, secondary P2 checks, ratio comparisons and scientific figures.
- `src/ion_noise.py`: default Stage 3 entry point, explicit Stage 1/2 options retained, configurable cache/grid/uncertainty options.
- `src/ion_params.py`: requested 0/100/200/400 microsecond readout sweep and uniform analysis grid; nominal hardware parameters unchanged.
- `src/steane_noise.py`: six-bin trigger histogram and no-trigger count in detailed MC reporting only; no RNG or quantum simulation changes.
- `src/ion_timing.py`: idle exposure per data qubit for all adaptive paths and their sampled average; existing timing/insertion logic unchanged.
- `tests/test_ion_calibration.py`: interpolation/threshold analytic checks, round trips, pooling, joint uncertainty, reproducibility, attribution conservation, and exposure accounting.
- `results/uniform_mc_calibration.json`: 11M uniform shots, with seeds, checks and pair comparisons/breakdowns.
- `results/measurement_duration_sweep.json`: four readout points, 3M shots each, their checks, trigger metrics and path exposures.
- `results/ion_noise_results.json`: all physical inputs, provenance, Monte Carlo, pair attribution, calibrated p_eff and ratios, confidence intervals, readout study, nominal exposure and marked Stage 2 T2 reference.
- `figures/ion_noise_breakdown.png`, `figures/ion_noise_calibration.png`: secondary type attribution and primary MC calibration, visually inspected.
- `docs/noise_model.md`, `README.md`: formulas, inference assumptions, uncertainty scope, cached runs and headline results.

Checks:
- Nominal Jones flagged single-fault failures: 0 / 1332; bare-only faults: 0 / 888, skipped on noiseless input; incoming errors: 0 / 42. Frame/tableau: 0 mismatches, including another 888 forced-trigger bare-fault agreement trials.
- Readout 0 microseconds: 0 / 1332 flagged, 0 / 882 bare-only, 0 / 42 incoming; 0 mismatches including forced bare execution. Three zero-length bare idle gaps are omitted (six fewer basis trials). Other readout points have nominal support and all checks pass.
- H2 timing-free checks remain 0 / 1152 flagged, 0 / 756 bare-only, 0 / 42 incoming; 0 simulator mismatches. Uniform channel support checks: 0 / 1128 flagged, 0 / 744 bare-only, 0 / 42 incoming; 0 mismatches.
- P1 = 0 for both presets, all readout points and uniform grid points.
- Full pytest: 43 passed in 134.82 seconds. Final fast suite: 38 passed in 36.84 seconds (5 slow deselected). Uniform and infinite-T2 exact regressions pass unchanged, including 542/3M uniform failures, 92840 malignant pairs, all A dictionaries, 0/1218 with flags, 48/1218 ignoring flags and 0/1176 baseline simulator mismatches.
- Confidence tests show that including both MC sources gives wider intervals than holding either source effectively fixed. Bootstrap is reproducible; extrapolation is rejected.
- Uniform interpolation round trip: maximum absolute error recorded below. The observed point grid is strictly increasing, with no point-estimate pooling.
- Nominal Stage 3 P1/P2 and all counting dictionaries exactly reproduce Stage 2; the fresh nominal MC differs by 0.908 combined standard errors.
- All secondary P2 p_eff values for the nominal presets lie inside the propagated primary MC 95% intervals. Raw P2 and survival-corrected P2 themselves remain attribution terms and need not equal p_L.
- Source/config fingerprints match all caches; each target has six explicit bracketing grid points at >=1M shots each. Artifact checks verify attribution sums and trigger/path accounting. No decoder, circuit gate order or adaptive stopping rule changed.
- `compileall` and `git diff --check` pass. Both figures were visually inspected. Default `python ion_noise.py` runs/resumes Stage 3 and writes the required final JSON and figures.

Calibration:
- 11 points: [0.00085, 0.00095, 0.00105, 0.00115, 0.00125, 0.00135, 0.00145, 0.00155, 0.00165, 0.00175, 0.00185], each 1,000,000 shots, distinct seeds 12026–12036.
- New readout experiments: four points × 3M shots, distinct seeds 3026/3126/3226/3326. Total new experimental sampling: 23M shots. H2 reuses the approved 3M-shot timing-free sample (seed 2026); no H2 timing is invented.
- Uniform MC pseudothreshold: 0.00133488898353; 95% CI [0.00129130961, 0.00145135925].
- Uniform round-trip max absolute error: 6.5052e-19.
- 20,000 joint binomial bootstrap draws, seed 302026; 1.725% of sampled grids required pooling; all 20,000 threshold/target draws remained bracketed. The point curve did not require pooling.

### Primary p_eff and ratios

All intervals below are 95% sampling intervals, propagated from target and uniform MC. Hardware and interpolation-model uncertainty are excluded.

| Preset | Primary MC p_L ± stderr | p_eff (MC) [95% CI] | p_eff / p_2q [95% CI] | p_eff / uniform pseudothreshold [95% CI] |
| --- | ---: | ---: | ---: | ---: |
| JONES_MURALI_2026 | 0.00204033333333 ± 2.6052322e-05 | 0.00167566902067 [0.001638817566, 0.001710165456] | 1.675669 [1.638817566, 1.710165456] | 1.2552872 [1.156312819, 1.306449499] |
| H2_LIKE (Stage 1 only) | 0.00087 ± 1.7021977e-05 | 0.0010734088134 [0.001032147376, 0.001135605607] | 0.83860064 [0.8063651373, 0.8871918802] | 0.80411841 [0.7380704684, 0.857423151] |

Jones–Murali lies above the one-gadget uniform pseudothreshold; timing-free H2 lies below it, with both ratio intervals remaining on their respective sides of one. This is a comparison under the defined uniform calibration, not a many-round threshold claim.

### Pair attribution and secondary check

| Preset | P1 | Raw P2 | Survival-corrected P2 | Raw / MC | Corrected / MC | P2-based p_eff | Within primary 95% CI? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| JONES_MURALI_2026 | 0 | 0.00218179149881 | 0.00182009035186 | 1.069331 | 0.8920554 | 0.0016855868979 | True |
| H2_LIKE | 0 | 0.00088465601772 | 0.000805656034421 | 1.016846 | 0.9260414 | 0.00107332622866 | True |

The JSON retains logical-X/Z, type-pair and stage-pair breakdowns for both raw and survival-corrected values. The figure splits mixed-type pair credit equally; credits sum to raw P2.

### Trigger metrics

| Preset | Trigger rate ± stderr | Share of failures through bare | Reset | Measurement | 2q | Idle | Mixed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| JONES_MURALI_2026 | 0.105150666667 ± 0.00017710073 | 6053/6121 = 0.9888907 | 0.05682066667 | 0.01135966667 | 0.02859266667 | 0.008126333333 | 0.0002513333333 |
| H2_LIKE | 0.0499503333333 ± 0.00012577135 | 2552/2610 = 0.97777778 | 0.0005803333333 | 0.01172433333 | 0.03761433333 | 0 | 3.133333333e-05 |

Exclusive rates sum to trigger rate. The JSON also includes overlapping contributing-type rates and equal-share mixed allocation, as defined in Stage 2. Every readout/grid point retains its trigger and bare-failure metrics.

### Jones–Murali readout-duration sensitivity

Only measurement duration changes. Measurement errors, MX rotation time, resets, gates and T2 remain nominal. Each point uses 3M shots.

| Measurement duration (μs) | MC p_L ± stderr | p_eff (MC) [95% CI] | Raw P2 | Survival-corrected P2 | Trigger rate | Failure share through bare |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.00167666666667 ± 2.3621004e-05 | 0.00153982523131 [0.001503013479, 0.001572012799] | 0.00181939901792 | 0.00153670991888 | 0.09957866667 | 0.9918489066 |
| 100 | 0.00179033333333 ± 2.4407158e-05 | 0.00158400704318 [0.001549512904, 0.001611361568] | 0.00190438991454 | 0.001603505754 | 0.1012036667 | 0.9912493018 |
| 200 | 0.00186066666667 ± 2.4881081e-05 | 0.00160894295754 [0.001580474676, 0.001640772799] | 0.00199311957308 | 0.00167301856298 | 0.1023916667 | 0.9888928699 |
| 400 | 0.00204033333333 ± 2.6052322e-05 | 0.00167566902067 [0.001638817566, 0.001710165456] | 0.00218179149881 | 0.00182009035186 | 0.1051506667 | 0.9888907041 |

400 μs − 0 μs logical-rate difference: 0.000363666666667 ± 3.5166394e-05; 95% CI [0.000295, 0.0004323333333]. Removing readout duration lowers nominal p_L by 17.8239%.

Relative to the Stage 1 no-idle reference p_L=0.001592, this difference is 81.1152% of the total nominal logical idle increment; 95% interval [68.1127%, 95.2460%]. This is a ratio of logical-rate increments with shared-numerator covariance propagated by resampling three experiments, not a microscopic failure allocation. All live-ion waiting changes, including ancillary waits; it is not solely a data-qubit intervention.

### Nominal Jones–Murali idle exposure per data qubit

The main values are actual-path MC averages over the fresh nominal 3M shots, rather than the sum over all potential circuits. Probability sums are expected idle-Z fault counts, not net Z probabilities. The ± values are sampling standard errors from the path histogram.

| Data qubit (0-indexed) | Idle exposure (seconds) ± stderr | Summed local Z probability ± stderr |
| ---: | ---: | ---: |
| 0 | 0.00691133266 ± 3.902102e-07 | 0.00157038308349 ± 8.866886e-08 |
| 1 | 0.00690994191333 ± 3.8929e-07 | 0.00157006754594 ± 8.846e-08 |
| 2 | 0.00682633998 ± 3.826338e-07 | 0.00155112393518 ± 8.694932e-08 |
| 3 | 0.00690848683333 ± 3.89651e-07 | 0.00156973741221 ± 8.854178e-08 |
| 4 | 0.0068248849 ± 3.830244e-07 | 0.00155079978794 ± 8.703806e-08 |
| 5 | 0.00682349415333 ± 3.821531e-07 | 0.00155049053227 ± 8.684084e-08 |
| 6 | 0.00673989222 ± 3.755823e-07 | 0.0015315675524 ± 8.535045e-08 |

Readout alone accounts for 0.00484848586667 seconds of waiting per data qubit on average (about 70–72% of their idle seconds). All seven conditional path exposure tables are in the JSON. For comparison, the no-trigger full flagged path has 0.006660–0.006820 seconds of data idle; the longest executed path has 0.010110–0.010430 seconds.

| Execution path | MC shots | Probability | Round duration (seconds) |
| --- | ---: | ---: | ---: |
| no_trigger | 2684548 | 0.8948493333 | 0.0069 |
| trigger_F0 | 48985 | 0.01632833333 | 0.00484 |
| trigger_F1 | 49422 | 0.016474 | 0.00599 |
| trigger_F2 | 48915 | 0.016305 | 0.00714 |
| trigger_F3 | 60146 | 0.02004866667 | 0.00829 |
| trigger_F4 | 54884 | 0.01829466667 | 0.00944 |
| trigger_F5 | 53100 | 0.0177 | 0.01059 |

Assumptions made:
- Uniform p_eff is a one-round MC equivalence under the same uniform channel conventions; it is not a hardware-averaged infidelity. H2 remains timing-free Stage 1 only. Jones gate rates and the H2 convention retain their earlier placeholder/unverified caveats.
- Confidence intervals propagate binomial target and calibration sampling, with common grid draws for threshold ratios. They are conditional on piecewise log-log interpolation and fixed physical input rates, excluding interpolation bias and hardware calibration uncertainty. Small sampled reversals are pooled; no extrapolation is used.
- The readout sweep is a hypothetical duration-only intervention: readout error stays finite even at zero duration, and MX basis rotations retain nominal duration/error. It changes data and ancillary waiting.
- Primary nominal Stage 3 Jones MC is a fresh independent 3M-shot sample, so its point estimate differs through sampling from the approved Stage 2 nominal sample. The Stage 2 T2 sweep remains a marked reference; this grid is not used to extrapolate p_eff for its high-error points.
- Summed per-gap idle-Z probabilities are expected physical fault counts, not a net probability after intervening gates. Exposure is averaged over actual executed paths; all potential locations are still used in the separate survival-corrected pair attribution, as before.
- Pair-type figures allocate half of a mixed-type pair to each type, preserving the raw P2 total. Monte Carlo remains primary.

Open questions for the human:
- No issue blocks this Stage 3 checkpoint. H2 timing remains deferred pending a sourced preset.
- Grid/target sampling intervals do not establish uncertainty in the hardware parameters or the interpolation model; those remain explicit modeling limitations.
- Stage 4 transport hooks have not started. Stop here for review.


## Stage 4 — Optional transport interface and synthetic validation

Scope: Stage 3 approved. Stage 4 only. The transport list remains explicit input;
all benchmark events are synthetic. No hardware schedule, trap layout, routing
policy or calibrated heating law is invented. The default CLI remains Stage 3.

Changes:
- `src/ion_transport.py`: event validation, liveness-aware waiting, next-executed-CX heating tables, branch-conditioned sampler/pair counter, explicit synthetic fixture and checkpoint runner.
- `src/ion_noise.py`: optional transport_events/heating_penalty inputs, context-aware fault support testing, Stage 4 CLI with empty default and explicit synthetic/custom event flags.
- `src/ion_params.py`: zero default penalty map and synthetic 0/0.0005/0.001 sensitivity values, marked placeholders. Sourced durations and other preset values unchanged.
- `src/steane_noise.py`: original independent count/MC bodies retained verbatim under private names; opt-in dispatch only for branch-dependent channels. Quantum simulators, decoder and circuits unchanged.
- `src/ion_timing.py`: exposure reporting accepts a supplied circuit map and extra durations; default timing/insertion behavior unchanged.
- `src/ion_cache.py`, `src/ion_analysis.py`: frozen semantic compatibility signature preserves approved Stage 3 data only when independent kernels, default dispatch, protocol/timing primitives, ordered baseline channels/circuits and study parameters match.
- `tests/test_ion_transport.py`: event kinds/liveness/indices, summed and consumed heating, cross-circuit/reset/adaptive-branch carry, exact conditional oracle, seeded sampler equivalence, empty regressions, synthetic FT and P2 monotonicity.
- `results/stage4_transport_results.json`: primary MC, attribution and comparisons, checks, explicit synthetic events, placeholder penalties, heated CX rates, trigger/failure metrics and exposure; Stage 3 baseline/H2 data retained.
- `README.md`, `docs/noise_model.md`: stable original-op indexing, branch semantics, default behavior, cache preservation and synthetic-only scope.

Checks:
- Empty transport: all circuits and ordered channels match Stage 3; P1/P2 and every counting dictionary match exactly. Replaying 3M shots at the approved seed 3326 exactly reproduces all Stage 3 MC metrics (6121 failures, 315452 triggers), including its trigger histogram.
- Rebuilding the Stage 3 final report using approved caches is byte-for-byte identical, including p_eff confidence intervals and ratios. No Stage 0–3 numerical results or figures were modified.
- Synthetic flagged single faults: 0 / 1344, including 192 idle-fault runs. Bare-only single faults: 0 / 912, including 156 idle-fault runs; these never execute on noiseless input and are reported separately. Incoming errors: 0 / 42.
- Synthetic frame/tableau cross-check: 0 mismatches / 2256 channel fault trials plus 42 incoming trials. Another 912 forced-trigger bare-fault trials per case agree, exercising actual bare execution; these are not added to the single-fault FT headline.
- The three synthetic heating levels have identical support and all pass the above checks; P1 = 0 at every level.
- Asymmetric branch-specific heating tested against an exact sparse oracle: next CX depends on the actual F stop; the MC and weighted pair counter reproduce the oracle. Support testing includes payloads available in any branch.
- Conditional and independent samplers give identical seeded metrics when branch tables coincide. Reset does not clear pending heating; accumulated penalties are consumed once at the next CX. Tests and explicit checks cover all four transport kinds, both CX participants, invalid inputs and probability overflow.
- P2 strictly increases with the synthetic heating penalty. This is deterministic pair attribution; finite-shot MC monotonicity is not imposed as a test.
- Full pytest: 61 passed in 267.51 s. Final fast suite: 55 passed in 49.84 s (6 slow deselected). Legacy uniform and infinite-T2 exact regressions all pass.
- `compileall`, `git diff --check`, artifact-sum and trigger/path accounting checks pass. H2 remains timing-free Stage 1 only.

Numbers: primary p_L is Monte Carlo; P2 is attribution. All new MC rows use 3M shots. The no-event row replays seed 3326; synthetic rows use independent seeds 6026/6027/6028. Heating increments are additive probability per event, not hardware estimates.

| Case | P1 | MC p_L ± stderr | Raw P2 | Survival-corrected P2 | Raw/MC | Corrected/MC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| No events (Stage 3 regression) | 0 | 0.00204033333333 ± 2.6052322e-05 | 0.00218179149881 | 0.00182009035186 | 1.069331 | 0.8920554 |
| Synthetic, shuttle penalty=0 | 0 | 0.00208233333333 ± 2.6318544e-05 | 0.00218859848822 | 0.00182533391532 | 1.051032 | 0.876581 |
| Synthetic, shuttle penalty=0.0005 | 0 | 0.00226233333333 ± 2.7430003e-05 | 0.00246016452556 | 0.00202844645241 | 1.087446 | 0.896617 |
| Synthetic, shuttle penalty=0.001 | 0 | 0.00256566666667 ± 2.9206643e-05 | 0.00274780167401 | 0.00224023112511 | 1.070989 | 0.8731575 |

The empty-input nominal p_eff and ratio confidence intervals are exactly those approved in Stage 3. H2 likewise keeps its approved timing-free 3M-shot primary MC and p_eff. Synthetic hook trials are not assigned a new hardware p_eff; no extrapolation of the Stage 3 grid is made.

### Trigger metrics

Rates below are probabilities per shot. Exclusive categories plus mixed sum to trigger rate; overlapping and equal-share splits are also retained in the JSON.

| Case | Trigger rate ± stderr | Failures through bare / total | Bare-failure share | Reset | Measurement | 2q | Idle | Mixed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| No events (Stage 3 regression) | 0.105150666667 ± 0.00017710073 | 6053 / 6121 | 0.9888907041 | 0.05682066667 | 0.01135966667 | 0.02859266667 | 0.008126333333 | 0.0002513333333 |
| Synthetic, shuttle penalty=0 | 0.105081 ± 0.00017704894 | 6185 / 6247 | 0.9900752361 | 0.05675166667 | 0.011327 | 0.02862033333 | 0.008131 | 0.000251 |
| Synthetic, shuttle penalty=0.0005 | 0.108053333333 ± 0.0001792371 | 6702 / 6787 | 0.9874760572 | 0.056628 | 0.01136566667 | 0.03170833333 | 0.008094 | 0.0002573333333 |
| Synthetic, shuttle penalty=0.001 | 0.111393333333 ± 0.00018164513 | 7597 / 7697 | 0.9870079252 | 0.05669933333 | 0.011305 | 0.03500366667 | 0.008098 | 0.0002873333333 |
| H2_LIKE, Stage 1 only | 0.0499503333333 ± 0.00012577135 | 2552 / 2610 | 0.9777777778 | 0.0005803333333 | 0.01172433333 | 0.03761433333 | 0 | 3.133333333e-05 |

Transport waiting is credited to idle and heated CX faults to 2q. No new causal category is introduced.

### Synthetic fixture and duration/heating bookkeeping

- Two A=7 shuttle markers per F/B circuit, before the first original CX and before syndrome readout: 24 potential events total. Each lasts the sourced 5 microseconds. This is an interface test, not routing or a claim of hardware feasibility.
- 270 potential fault locations: 96 original operations and 174 idle locations (96 flagged, 78 bare). New ancilla gaps created by transport are included in both simulators.
- Each executed circuit gains 10 microseconds. No-trigger F0..F5 gains 60 microseconds; the longest executed F0..F5+B0..B5 path gains 120 microseconds. Skipped events do not contribute to actual timing or realized heating.
- At penalty h, F0 first CX gets p_2q+h; later first CXs get p_2q+2h because the preceding out-marker and current in-marker accumulate. All other CXs retain nominal p_2q. The final out-marker has no later gate within this round.
- This symmetric fixture has equal B rate tables across trigger branches and therefore uses the unchanged independent kernel. Asymmetric fixtures in tests exercise the conditional implementation and verify correct carry into B versus a later F circuit.

### Attribution — No events (Stage 3 regression)

| Location type pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| 2q*2q | 0.000388844444444 | 0.000323602954618 |
| 2q*idle | 0.000244229987957 | 0.00020308726939 |
| 2q*meas | 0.000178139673778 | 0.000148251372267 |
| 2q*reset | 0.000888268146667 | 0.0007422074391 |
| idle*idle | 6.2034119455e-05 | 5.15425354296e-05 |
| idle*meas | 3.40953488086e-05 | 2.83514271743e-05 |
| idle*reset | 0.000169702130754 | 0.000141680701929 |
| meas*meas | 6.04005280053e-06 | 5.02666828385e-06 |
| meas*reset | 6.0239463472e-05 | 5.03342626347e-05 |
| reset*reset | 0.00015019813068 | 0.00012600572103 |

| Stage pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| FB | 0.0019670291389 | 0.00164135880758 |
| FF | 0.000214762359919 | 0.000178731544274 |

### Attribution — Synthetic, shuttle penalty=0

| Location type pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| 2q*2q | 0.000388844444444 | 0.000323527955278 |
| 2q*idle | 0.000247149479354 | 0.000205467533899 |
| 2q*meas | 0.000178139673778 | 0.000148217013016 |
| 2q*reset | 0.000888268146667 | 0.000742035422536 |
| idle*idle | 6.34565433849e-05 | 5.2712290623e-05 |
| idle*meas | 3.45076536685e-05 | 2.86876520035e-05 |
| idle*reset | 0.000171754899977 | 0.000143361430175 |
| meas*meas | 6.04005280053e-06 | 5.02550328582e-06 |
| meas*reset | 6.0239463472e-05 | 5.03225969919e-05 |
| reset*reset | 0.00015019813068 | 0.000125976517509 |

| Stage pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| FB | 0.00197246542093 | 0.00164550420796 |
| FF | 0.000216133067296 | 0.000179829707363 |

### Attribution — Synthetic, shuttle penalty=0.0005

| Location type pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| 2q*2q | 0.00052899111111 | 0.000435334195601 |
| 2q*idle | 0.000295266357579 | 0.000242737069885 |
| 2q*meas | 0.000192057792889 | 0.000157989969675 |
| 2q*reset | 0.00095765252 | 0.000790949351035 |
| idle*idle | 6.34565433849e-05 | 5.21086778082e-05 |
| idle*meas | 3.45076536685e-05 | 2.83591473194e-05 |
| idle*reset | 0.000171754899977 | 0.00014171978654 |
| meas*meas | 6.04005280053e-06 | 4.9679558306e-06 |
| meas*reset | 6.0239463472e-05 | 4.97463487573e-05 |
| reset*reset | 0.00015019813068 | 0.000124533949951 |

| Stage pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| FB | 0.00218726194229 | 0.00180387865681 |
| FF | 0.000272902583274 | 0.000224567795593 |

### Attribution — Synthetic, shuttle penalty=0.001

| Location type pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| 2q*2q | 0.000685208888888 | 0.000557910666916 |
| 2q*idle | 0.000343383235804 | 0.000279199907708 |
| 2q*meas | 0.000205975912 | 0.00016753975679 |
| 2q*reset | 0.00102703689333 | 0.000838746285138 |
| idle*idle | 6.34565433849e-05 | 5.15113951951e-05 |
| idle*meas | 3.45076536685e-05 | 2.80340877261e-05 |
| idle*reset | 0.000171754899977 | 0.000140095359132 |
| meas*meas | 6.04005280053e-06 | 4.91101188642e-06 |
| meas*reset | 6.0239463472e-05 | 4.91761437467e-05 |
| reset*reset | 0.00015019813068 | 0.000123106510872 |

| Stage pair | Raw P2 | Survival-corrected contribution |
| --- | ---: | ---: |
| FB | 0.00241386290809 | 0.00196839420282 |
| FF | 0.000333938765919 | 0.000271836922287 |

### Data idle exposure in the time-only synthetic fixture

MC-weighted actual-path means; the JSON also contains all seven conditional path tables and both nonzero-heating cases. Probability sums are expected Z-fault counts, not net Z probabilities.

| Data qubit | Idle seconds | Summed Z probability |
| ---: | ---: | ---: |
| 0 | 0.00697597249667 | 0.00158506787318 |
| 1 | 0.00697457499 | 0.00158475082346 |
| 2 | 0.00689096785667 | 0.00156580727029 |
| 3 | 0.00697312624333 | 0.00158442214901 |
| 4 | 0.00688951911 | 0.00156548458182 |
| 5 | 0.00688812160333 | 0.0015651738162 |
| 6 | 0.00680451447 | 0.001546250896 |

Assumptions made:
- Indices are original stabilizer ops before idle insertion; the terminal index is allowed. Events at a gap execute sequentially. Only live ions receive dephasing; an unprepared/measured moved ion still accumulates deterministic heating for its future CX.
- Movement is logically ideal apart from added idle Z and the configurable next-CX depolarizing increment. Penalties sum over events and over both participants, survive resets, and are consumed once. Unused increments at the end of this one-round model have no realized gate effect; no cross-round thermal state is modeled.
- Bare heating rates follow actual executed flagged prefixes. Latent unexecuted F/B faults are still sampled for consistency with the full-pool survival attribution; they do not affect actual execution. Exact survival products use the branch's rates; a branch-specific mean approximation is also retained.
- All nonzero heating values are synthetic placeholders, not calibrated heating physics. No gate duration, decoder, gate order, adaptive rule, native synthesis or hardware transport policy changes.
- MC remains primary. P2 raw/survival contributions, trigger causes, failure-through-bare share and stage-separated single-fault counts are reported as requested. Nominal p_eff is unchanged; synthetic test cases are not used to make a hardware p_eff claim.
- Empty-input semantic cache compatibility is explicit and validated; approved Stage 0–3 artifacts remain unchanged. The user's untracked `STAGE_LOG0.md` was left untouched.

Open questions for the human:
- None blocks this Stage 4 checkpoint. Actual schedules, transport locations/feasibility and nonzero heating penalties still require external scheduler/hardware inputs.
- H2 transport/idle modeling remains deferred until a sourced timing preset is supplied.
- Stop here. No hardware schedule or work beyond Stage 4 has started.
