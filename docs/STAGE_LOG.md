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


## Stage A1 — Parameterized flag check order and timing policy

Changes:
- `src/flag_design.py`: optional validated `round_order` for both rounds, original kind/support decoder selection, SEQ/PAR_MR timing and live-qubit idle insertion, actual-path exposure, and unchanged-kernel adapters.
- `src/ion_transport.py`: optional circuit map and idle inserter; default Stage 4 path unchanged. Event positions follow the executed order, and event insertions split parallel layers.
- `src/ion_params.py`: centralized user-requested design scan grids; no new hardware values.
- `src/check_order.py`: A1-only correctness matrix and cached baseline-order timing references; atomic result saves and source/configuration checks.
- `tests/test_flag_design.py`, `tests/fixtures/stages0_4_artifacts.json`: order, timing, conditional synthetic transport, default replay, exhaustive-matrix and frozen-artifact checks.
- `README.md`, `docs/noise_model.md`: usage, timing assumption, transport-index convention, and stage scope.
- `results/stageA1_fault_checks.json`, `stageA1_reference_SEQ.json`, `stageA1_reference_PAR_MR.json`, `stageA1_validation.json`: new checkpoint outputs. No A2/B output generated.

Frozen-result checks:
- All 16 approved `results/` and `figures/` files are SHA-256 identical to commit `d298328beb86d96699a2f3c8804eff3268419e9a`; manifest records each digest.
- Default-order functions delegate to the approved instrument; original circuit and decoder source files are unchanged. Approved Stage 3 semantic fingerprint remains `481c4d8668af7a000bfed3b2e92163228375fd8a36d4af386dea49d5e8539dec`.
- Full suite: **81 passed in 401.31 s** (61 existing + 20 new, including all slow checks); `git diff --check` passed.
- Flagged FT regression: 0/1218; ignoring flags: 48/1218; original frame/tableau regression: 0/1176 — all pass.
- Uniform malignant-pair total 92840 and all A dictionaries — exact pass; uniform seed 2026 at 3M shots: exactly 542 failures.
- Default-order Jones seed 3326 at 3M shots: exactly 6121 failures, with the full MC/trigger dictionary equal to approved Stage 3 and Stage 4 empty transport; default pair counts equal Stage 3.
- Stage 4 empty-transport channels, sampling, pair counts, cached Stage 3 report reconstruction and synthetic P2 monotonicity — pass. No approved artifact was rewritten.

FT checks: all four orders (`ZZZXXX`, `XXXZZZ`, `ZXZXZX`, `XZXZXZ`) pass every row below. Counts are **per order and per scan value**, not pooled. All failure numerators are zero. Positive-support-identical scans share exhaustive results: 36 reported rows, 20 distinct exhaustive geometries. Bare-only faults are tested separately from flagged faults; they never execute on noiseless input. Incoming tests cover 7 qubits × X/Y/Z × both logical bases.

| Preset | Policy | Flagged failures/runs | Bare-only failures/runs | Incoming failures/runs | Frame/tableau mismatches/runs | Additional forced-trigger bare checks |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| UniformNoise p={0.0005,0.001,0.0015} | timing-free | 0/1128 | 0/744 | 0/42 | 0/1872 | 0/744 |
| TWO_QUBIT_ONLY p={0.0005,0.001,0.002} | timing-free | 0/1080 | 0/720 | 0/42 | 0/1800 | 0/720 |
| H2_LIKE (Stage 1 only) | timing-free | 0/1152 | 0/756 | 0/42 | 0/1908 | 0/756 |
| JONES_MURALI_2026 | SEQ | 0/1332 | 0/888 | 0/42 | 0/2220 | 0/888 |
| JONES_MURALI_2026 | PAR_MR | 0/1326 | 0/888 | 0/42 | 0/2214 | 0/888 |

Cross-check: zero mismatches throughout. Tests also reverse support-row order, compare incoming-error syndromes/corrections between simulators, preserve noiseless states in both bases, and exercise an asymmetric late synthetic transport event with prefix-dependent heating and forced bare-round faults. Transport is still synthetic.

Numbers: **baseline order only**, nominal Jones–Murali, empty transport. These are the required SEQ/PAR_MR timing reference rows for the same gadget, not an A2 order-design comparison. Monte Carlo is primary, 5,000,000 shots per arm (independent seeds 73001/73002); P1=0 for both.

| Policy | MC failures | MC p_L ± stderr | Raw P2 | Exact survival-corrected P2 | Raw/MC | Corrected/MC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| SEQ | 10200 | 0.00204 ± 2.01783964e-05 | 0.00218179149881 | 0.00182009035186 | 1.069506 | 0.892201 |
| PAR_MR | 9551 | 0.0019102 ± 1.95271664e-05 | 0.00201683134625 | 0.00169104430217 | 1.055822 | 0.885271 |

PAR_MR minus SEQ: Δp_L = -0.0001298 ± 0.00002807985 (combined stderr), -4.6225σ; significant for this timing assumption at the specified ≥3σ rule. Combined stderr is 1.4700% of the smaller p_L, meeting the ≤2% criterion. ΔP2 = -0.000164960153. Raw P2 is above MC by 6.95%/5.58% (SEQ/PAR_MR); the survival-corrected estimate is below by 10.78%/11.47%. Neither replaces adaptive MC. The exact correction multiplies each pair weight by the product over all other potential locations; JSON also retains the mean-probability approximation and type/stage attribution.

Trigger rates are probabilities per shot. The exclusive causes below sum to the total; mixed means multiple contributing types, with fractional and overlapping contributor rates also saved in JSON. Failure-through-bare is conditional on an MC failure.

| Policy | Total trigger | Reset | Measurement | 2q | Idle | Mixed | Failures through bare / all failures |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| SEQ | 0.1049754 | 0.0568204 | 0.011358 | 0.028475 | 0.0080622 | 0.0002598 | 10106/10200 (99.07843%) |
| PAR_MR | 0.101814 | 0.0570394 | 0.0112968 | 0.0285042 | 0.0047432 | 0.0002304 | 9454/9551 (98.98440%) |

Per-position trigger rates (F0…F5, all-shot denominator):

- SEQ: 0.0162038, 0.0163502, 0.0162802, 0.0199438, 0.0184906, 0.0177068.
- PAR_MR: 0.0161858, 0.0163438, 0.0161172, 0.0187716, 0.0174884, 0.0169072.

Timing: each flagged check is 1150 us SEQ / 700 us PAR_MR; all-six-flagged no-trigger round is 6.9 ms / 4.2 ms. Bare Z/X checks remain 610/620 us under either policy. MC-weighted per-data-qubit readout idle is 0.00484903008 s SEQ / 0.002573230203 s PAR_MR. Full conditional paths, round durations and exposure errors are in the JSON.

| Data qubit | SEQ mean idle s | SEQ summed Z probability | PAR_MR mean idle s | PAR_MR summed Z probability |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 0.006912069708 | 0.00157055053742 | 0.004317464816 | 0.000981097382033 |
| 1 | 0.006910676076 | 0.00157023434523 | 0.004316111528 | 0.00098079002791 |
| 2 | 0.006827070948 | 0.00155129001833 | 0.004232663864 | 0.000961852964404 |
| 3 | 0.006909230172 | 0.00156990629339 | 0.004314713232 | 0.000980472451725 |
| 4 | 0.006825625044 | 0.00155096795216 | 0.004231265568 | 0.000961538540756 |
| 5 | 0.006824231412 | 0.00155065804283 | 0.00422991228 | 0.000961234483198 |
| 6 | 0.006740626284 | 0.00153173434939 | 0.004146464616 | 0.000942304242141 |

Assumptions made:
- PAR_MR groups adjacent disjoint resets or measurements, including the existing composite RX/MX durations. This is a timing assumption, not a demonstrated device schedule; gates are sequential. It changes only waiting/idle noise. No hardware values or placeholders were added.
- Bare syndrome rows retain their original support addresses, even if support order is also permuted. The triggered flagged lookup uses actual kind and support; decoders and adaptive rules are unchanged.
- Transport indices identify positions in the executed order. Insertions are serial barriers; branch-dependent heating uses only the executed flagged prefix. No real routing or layout is introduced.
- Uniform and two-qubit-only scan probabilities are analysis inputs from the prompt. H2 remains timing-free, Stage 1 only. Jones uses nominal sourced timing and empty transport.
- Monte Carlo remains primary; P2 and survival-corrected P2 are attribution. Baseline timing references are included here as requested. Order MC, logical-X/Z comparisons, p_eff and calibration extension belong to A2 and have not started.

Open questions for the human:
- None blocks A1. Review the order/timing interfaces and PAR_MR assumption before A2.
- Stop here. Stage A2 and all Stage B work await approval.


## Stage A2 — Stabilizer-check order experiment

Changes:
- `src/order_experiment.py`: 36-row experiment, exact pair attribution/reweighting, independent-stream MC allocation, resumable atomic caches, extended calibration, same-policy comparisons and nine-panel figure.
- `src/design_mc.py`: exact sparse batched independent-location Monte Carlo, original decoder tables, four linear trigger-attribution component frames and logical-X/Z split. No approved scalar kernel changes.
- `src/ion_params.py`: seven lower analysis-only uniform calibration points, no hardware parameters.
- `src/check_order.py`: A1 fingerprint compatibility for only the appended literal A2 grid assignment; runtime overrides are never excluded.
- `tests/test_design_mc.py`, `tests/test_order_experiment.py`, `tests/fixtures/stageA1_artifacts.json`: payload/mixed-fault scalar/tableau oracles, sampling probabilities, resumability, pair reweighting and four-file A1 freeze manifest.
- `results/stageA_check_order.json`, `results/stageA2_validation.json`, `results/stageA2_cache/`, `figures/check_order_comparison.png`: complete results, validation, reproducible caches and visually inspected figure.
- `README.md`, `docs/noise_model.md`, this log: instructions, statistical methods, assumptions and results.

Frozen-result checks:
- All 16 approved Stage 0–4 result/figure files and all four approved A1 JSON files remain SHA-256 byte-identical. Original flag circuit, decoder, scalar MC and timing builders are unchanged.
- Validation covers all **99 distinct test nodes**: final fast run **91 passed, 8 slow deselected in 52.95 s**; all **eight slow nodes passed in the 95-test full run (427.58 s)**. Focused A2 run: 18 passed. `git diff --check` passes.
- Original flagged FT: 0/1218; ignoring flags: 48/1218; original frame/tableau: 0/1176 — pass.
- Uniform 92840 malignant pairs and A dictionaries — exact pass; uniform seed 2026: exactly 542/3M failures.
- Jones seed 3326: exactly 6121/3M failures; full MC/trigger dictionary and default pair counts match Stage 3/empty Stage 4. Stage 4 empty replay, report reconstruction and synthetic monotonicity regressions — pass.
- Independent new baseline MC agrees with approved A1: +0.7365σ under SEQ and -0.7152σ under PAR_MR, using combined MC errors. Frozen scalar RNG streams remain intact.

FT checks: circuit/channel supports are identical to A1. Every P1 and single-fault failure count remains zero. The approved 36-row A1 matrix is checked against its unchanged source signature and artifact hashes. Counts below apply to every one of the four orders and each scan value; bare-only faults are separate because they never execute on noiseless input.

| Preset | Policy | Flagged failures/runs | Bare-only failures/runs | Incoming failures/runs | Frame/tableau mismatches/runs | Forced-trigger bare cross-check |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Uniform p={0.0005,0.001,0.0015} | timing-free | 0/1128 | 0/744 | 0/42 | 0/1872 | 0/744 |
| TWO_QUBIT_ONLY p={0.0005,0.001,0.002} | timing-free | 0/1080 | 0/720 | 0/42 | 0/1800 | 0/720 |
| H2_LIKE, Stage 1 only | timing-free | 0/1152 | 0/756 | 0/42 | 0/1908 | 0/756 |
| JONES_MURALI_2026 | SEQ | 0/1332 | 0/888 | 0/42 | 0/2220 | 0/888 |
| JONES_MURALI_2026 | PAR_MR | 0/1326 | 0/888 | 0/42 | 0/2214 | 0/888 |

Cross-check: the new batch engine replays every positive single-fault payload and forced-trigger bare payload against the original scalar simulator across all 20 distinct support geometries. Mixed-fault tests include cancellations, all categories and changed flagged corrections; 50 explicit multi-fault shots per order/policy are also checked with tableau in both bases (800 both-basis trials). Zero frame/log/trigger/attribution/residual/logical-decision mismatches in the completed validation. Sparse sampling tests check independent locations, uniform and unequal payload probabilities, and deterministic cache resumption.

Numbers: 932,000,000 new target MC shots across 36 rows, plus 35,000,000 new calibration shots. All four Jones orders run under **both SEQ and PAR_MR**, nominal sourced parameters and empty transport. The other seven families remain timing-free; H2 remains **Stage 1 only**. No A2 comparison crosses timing policies.

### Findings

- For Jones–Murali, **XXXZZZ improves p_L under both timing policies**. It is the only tested order with a significant improvement over the same-policy baseline. SEQ improves by 9.94%; PAR_MR by 7.07% (point reductions). Exact differences and errors are below.
- Both alternating orders **increase** p_L in all nine families, each at >3σ. Raw P2 differences have the same sign. Interleaving does not improve this unchanged adaptive gadget/decoder.
- Under all uniform/CX-only rates, XXXZZZ versus ZZZXXX is a null result (|σ|≤0.514). H2's XXXZZZ difference is -1.782σ, also below the 3σ rule. No improvement claim is made for these seven families.
- The first-X trigger spike is reduced by X-first ordering: in SEQ the baseline F3 rate is 0.0200145, while XXXZZZ's X positions F0–F2 are 0.0163193, 0.0173970, 0.0176290. Under PAR_MR, baseline F3 is 0.0187693, while X-first F0–F2 are 0.0163158, 0.0169565, 0.0170650. Histograms below retain every position. Earlier X detection alone does not predict p_L: alternating orders also move X checks earlier but perform worse.
- The logical-Z-only share among Jones MC failures decreases from 58.55% to 51.43% under SEQ and 57.28% to 49.81% under PAR_MR when moving to XXXZZZ. Both-XZ failures are reported separately, not double-counted. These are descriptive shares, not additional significance claims.

| Jones policy | Baseline MC p_L ± stderr | XXXZZZ MC p_L ± stderr | Δp_L ± combined stderr | σ | ΔP2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| SEQ | 0.002060166667 ± 1.85109e-05 | 0.001855333333 ± 1.75684e-05 | -0.000204833333 ± 2.55206e-05 | -8.0262 | -0.0002187661444 |
| PAR_MR | 0.001891333333 ± 1.77377e-05 | 0.001757666667 ± 1.71006e-05 | -0.000133666667 ± 2.46385e-05 | -5.4251 | -0.0001484901804 |

### All MC comparisons

Differences are order minus ZZZXXX in the same preset/policy. Error is combined **one-sigma** stderr, not the 95% bars used in the figure. Significance is |σ|≥3 per comparison, without a multiple-comparison/family-wise claim. Largest combined stderr / smaller p_L = 1.4719%, below 2%. Shot counts were fixed initially using deterministic P2 with a variance margin; no MC-driven extension was needed.

| Preset / policy / order | M shots | MC p_L ± stderr | Δp_L ± combined stderr | σ | ≥3σ |
| --- | ---: | ---: | ---: | ---: | --- |
| UniformNoise(0.0005) / timing-free / ZZZXXX | 53 | 0.000188754717 ± 1.887e-06 | 0 ± 0 | 0.000 | no |
| UniformNoise(0.0005) / timing-free / XXXZZZ | 53 | 0.000189113208 ± 1.8888e-06 | 3.584906e-07 ± 2.6699e-06 | 0.134 | no |
| UniformNoise(0.0005) / timing-free / ZXZXZX | 53 | 0.000217075472 ± 2.0236e-06 | 2.832075e-05 ± 2.7669e-06 | 10.236 | yes |
| UniformNoise(0.0005) / timing-free / XZXZXZ | 53 | 0.000214245283 ± 2.0103e-06 | 2.549057e-05 ± 2.7572e-06 | 9.245 | yes |
| UniformNoise(0.001) / timing-free / ZZZXXX | 14 | 0.000742 ± 7.2774e-06 | 0 ± 0 | 0.000 | no |
| UniformNoise(0.001) / timing-free / XXXZZZ | 14 | 0.000743214286 ± 7.2834e-06 | 1.214286e-06 ± 1.0296e-05 | 0.118 | no |
| UniformNoise(0.001) / timing-free / ZXZXZX | 14 | 0.000854428571 ± 7.8089e-06 | 0.0001124286 ± 1.0674e-05 | 10.533 | yes |
| UniformNoise(0.001) / timing-free / XZXZXZ | 14 | 0.000851285714 ± 7.7945e-06 | 0.0001092857 ± 1.0664e-05 | 10.248 | yes |
| UniformNoise(0.0015) / timing-free / ZZZXXX | 6 | 0.00164983333 ± 1.6569e-05 | 0 ± 0 | 0.000 | no |
| UniformNoise(0.0015) / timing-free / XXXZZZ | 6 | 0.00165083333 ± 1.6574e-05 | 1e-06 ± 2.3435e-05 | 0.043 | no |
| UniformNoise(0.0015) / timing-free / ZXZXZX | 6 | 0.0018945 ± 1.7753e-05 | 0.0002446667 ± 2.4283e-05 | 10.076 | yes |
| UniformNoise(0.0015) / timing-free / XZXZXZ | 6 | 0.00186683333 ± 1.7623e-05 | 0.000217 ± 2.4188e-05 | 8.971 | yes |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZZZXXX | 103 | 9.60873786e-05 ± 9.6581e-07 | 0 ± 0 | 0.000 | no |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XXXZZZ | 103 | 9.64271845e-05 ± 9.6752e-07 | 3.398058e-07 ± 1.3671e-06 | 0.249 | no |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZXZXZX | 103 | 0.000105533981 ± 1.0122e-06 | 9.446602e-06 ± 1.399e-06 | 6.752 | yes |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XZXZXZ | 103 | 0.000106058252 ± 1.0147e-06 | 9.970874e-06 ± 1.4008e-06 | 7.118 | yes |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZZZXXX | 26 | 0.000381923077 ± 3.8319e-06 | 0 ± 0 | 0.000 | no |
| TWO_QUBIT_ONLY(0.001) / timing-free / XXXZZZ | 26 | 0.000380615385 ± 3.8254e-06 | -1.307692e-06 ± 5.4145e-06 | -0.242 | no |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZXZXZX | 26 | 0.000413961538 ± 3.9894e-06 | 3.203846e-05 ± 5.5316e-06 | 5.792 | yes |
| TWO_QUBIT_ONLY(0.001) / timing-free / XZXZXZ | 26 | 0.000421192308 ± 4.024e-06 | 3.926923e-05 ± 5.5567e-06 | 7.067 | yes |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZZZXXX | 7 | 0.00148285714 ± 1.4544e-05 | 0 ± 0 | 0.000 | no |
| TWO_QUBIT_ONLY(0.002) / timing-free / XXXZZZ | 7 | 0.00149342857 ± 1.4595e-05 | 1.057143e-05 ± 2.0605e-05 | 0.513 | no |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZXZXZX | 7 | 0.00159857143 ± 1.51e-05 | 0.0001157143 ± 2.0965e-05 | 5.519 | yes |
| TWO_QUBIT_ONLY(0.002) / timing-free / XZXZXZ | 7 | 0.00160742857 ± 1.5141e-05 | 0.0001245714 ± 2.0995e-05 | 5.933 | yes |
| H2_LIKE / timing-free / ZZZXXX | 12 | 0.000864416667 ± 8.4837e-06 | 0 ± 0 | 0.000 | no |
| H2_LIKE / timing-free / XXXZZZ | 12 | 0.000843166667 ± 8.3788e-06 | -2.125e-05 ± 1.1924e-05 | -1.782 | no |
| H2_LIKE / timing-free / ZXZXZX | 12 | 0.000963833333 ± 8.9578e-06 | 9.941667e-05 ± 1.2338e-05 | 8.058 | yes |
| H2_LIKE / timing-free / XZXZXZ | 12 | 0.000960166667 ± 8.9408e-06 | 9.575e-05 ± 1.2325e-05 | 7.769 | yes |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 6 | 0.00206016667 ± 1.8511e-05 | 0 ± 0 | 0.000 | no |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 6 | 0.00185533333 ± 1.7568e-05 | -0.0002048333 ± 2.5521e-05 | -8.026 | yes |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 6 | 0.00229 ± 1.9514e-05 | 0.0002298333 ± 2.6897e-05 | 8.545 | yes |
| JONES_MURALI_2026 / SEQ / XZXZXZ | 6 | 0.0022375 ± 1.9289e-05 | 0.0001773333 ± 2.6735e-05 | 6.633 | yes |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 6 | 0.00189133333 ± 1.7738e-05 | 0 ± 0 | 0.000 | no |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 6 | 0.00175766667 ± 1.7101e-05 | -0.0001336667 ± 2.4638e-05 | -5.425 | yes |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 6 | 0.00214316667 ± 1.8879e-05 | 0.0002518333 ± 2.5905e-05 | 9.722 | yes |
| JONES_MURALI_2026 / PAR_MR / XZXZXZ | 6 | 0.00208166667 ± 1.8607e-05 | 0.0001903333 ± 2.5707e-05 | 7.404 | yes |

### Secondary pair attribution, including survival correction

P1=0 for all rows. Raw P2 and exact survival-corrected P2 remain attribution; neither is used as the primary logical error estimate. Both are compared to MC below; discrepancies in MC stderr and the mean-probability survival approximation are also in JSON. Survival is the product over all other potential locations, including latent unexecuted locations, matching the approved convention.

| Preset / policy / order | Raw P2 | Corrected P2 | ΔP2 vs baseline | Δcorrected P2 | Raw/MC | Corrected/MC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| UniformNoise(0.0005) / timing-free / ZZZXXX | 0.0001919777778 | 0.0001831614255 | 0 | 0 | 1.01708 | 0.97037 |
| UniformNoise(0.0005) / timing-free / XXXZZZ | 0.0001919777778 | 0.0001831614255 | 0 | 0 | 1.01515 | 0.96853 |
| UniformNoise(0.0005) / timing-free / ZXZXZX | 0.0002190366667 | 0.000208977667 | 2.705889e-05 | 2.581624e-05 | 1.00903 | 0.96270 |
| UniformNoise(0.0005) / timing-free / XZXZXZ | 0.0002190366667 | 0.000208977667 | 2.705889e-05 | 2.581624e-05 | 1.02236 | 0.97541 |
| UniformNoise(0.001) / timing-free / ZZZXXX | 0.0007679111111 | 0.0006989833725 | 0 | 0 | 1.03492 | 0.94203 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 0.0007679111111 | 0.0006989833725 | 0 | 0 | 1.03323 | 0.94049 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 0.0008761466667 | 0.0007975036993 | 0.0001082356 | 9.852033e-05 | 1.02542 | 0.93338 |
| UniformNoise(0.001) / timing-free / XZXZXZ | 0.0008761466667 | 0.0007975036993 | 0.0001082356 | 9.852033e-05 | 1.02920 | 0.93682 |
| UniformNoise(0.0015) / timing-free / ZZZXXX | 0.0017278 | 0.001500416996 | 0 | 0 | 1.04726 | 0.90944 |
| UniformNoise(0.0015) / timing-free / XXXZZZ | 0.0017278 | 0.001500416996 | 0 | 0 | 1.04662 | 0.90888 |
| UniformNoise(0.0015) / timing-free / ZXZXZX | 0.00197133 | 0.001711897811 | 0.00024353 | 0.0002114808 | 1.04055 | 0.90361 |
| UniformNoise(0.0015) / timing-free / XZXZXZ | 0.00197133 | 0.001711897811 | 0.00024353 | 0.0002114808 | 1.05598 | 0.91701 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZZZXXX | 9.721111111e-05 | 9.4431789e-05 | 0 | 0 | 1.01169 | 0.98277 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XXXZZZ | 9.721111111e-05 | 9.4431789e-05 | 0 | 0 | 1.00813 | 0.97931 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZXZXZX | 0.00010587 | 0.0001028431152 | 8.658889e-06 | 8.411326e-06 | 1.00318 | 0.97450 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XZXZXZ | 0.00010587 | 0.0001028431152 | 8.658889e-06 | 8.411326e-06 | 0.99823 | 0.96969 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZZZXXX | 0.0003888444444 | 0.0003669223916 | 0 | 0 | 1.01812 | 0.96072 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XXXZZZ | 0.0003888444444 | 0.0003669223916 | 0 | 0 | 1.02162 | 0.96402 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZXZXZX | 0.00042348 | 0.0003996052833 | 3.463556e-05 | 3.268289e-05 | 1.02299 | 0.96532 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XZXZXZ | 0.00042348 | 0.0003996052833 | 3.463556e-05 | 3.268289e-05 | 1.00543 | 0.94875 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZZZXXX | 0.001555377778 | 0.001384864506 | 0 | 0 | 1.04891 | 0.93392 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XXXZZZ | 0.001555377778 | 0.001384864506 | 0 | 0 | 1.04148 | 0.92731 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZXZXZX | 0.00169392 | 0.001508218593 | 0.0001385422 | 0.0001233541 | 1.05965 | 0.94348 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XZXZXZ | 0.00169392 | 0.001508218593 | 0.0001385422 | 0.0001233541 | 1.05381 | 0.93828 |
| H2_LIKE / timing-free / ZZZXXX | 0.0008846560177 | 0.0008056560344 | 0 | 0 | 1.02341 | 0.93202 |
| H2_LIKE / timing-free / XXXZZZ | 0.0008844219929 | 0.0008054429328 | -2.340248e-07 | -2.131016e-07 | 1.04893 | 0.95526 |
| H2_LIKE / timing-free / ZXZXZX | 0.0009913073944 | 0.0009027778934 | 0.0001066514 | 9.712186e-05 | 1.02850 | 0.93665 |
| H2_LIKE / timing-free / XZXZXZ | 0.0009910601851 | 0.0009025527861 | 0.0001064042 | 9.689675e-05 | 1.03218 | 0.94000 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 0.002181791499 | 0.001820090352 | 0 | 0 | 1.05904 | 0.88347 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 0.001963025354 | 0.001637958409 | -0.0002187661 | -0.0001821319 | 1.05804 | 0.88284 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 0.002445684676 | 0.002040416984 | 0.0002638932 | 0.0002203266 | 1.06798 | 0.89101 |
| JONES_MURALI_2026 / SEQ / XZXZXZ | 0.002371951888 | 0.001979046748 | 0.0001901604 | 0.0001589564 | 1.06009 | 0.88449 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 0.002016831346 | 0.001691044302 | 0 | 0 | 1.06635 | 0.89410 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 0.001868341166 | 0.001566829377 | -0.0001484902 | -0.0001242149 | 1.06297 | 0.89143 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 0.002295578799 | 0.001924882382 | 0.0002787475 | 0.0002338381 | 1.07112 | 0.89815 |
| JONES_MURALI_2026 / PAR_MR / XZXZXZ | 0.0022461323 | 0.001883528673 | 0.000229301 | 0.0001924844 | 1.07901 | 0.90482 |

### MC effective rates and uniform pseudothreshold

The original 11-point, ≥1M-shot Stage 3 calibration is unchanged. Seven lower points at p={0.00015,0.00025,0.00035,0.00045,0.00055,0.00065,0.00075} have 5M shots each. All targets and their bootstrap intervals are bracketed, with no extrapolation. Each reported neighborhood has at least six points with ≥1M shots. The primary p_eff uses log-log MC interpolation and a 20,000-draw joint target/calibration binomial bootstrap. Grid draws are shared for threshold ratios. Individual ratio CIs and bracketing neighborhoods are in JSON.

Uniform MC pseudothreshold p_L=p: **0.00133488898353**, 95% CI **[0.00129094142128, 0.00145101331547]**. P2-based p_eff uses the original uniform A only as a secondary check.

| Preset / policy / order | MC p_eff [95% CI] | p_eff/p_2q | p_eff/pseudothreshold | P2 p_eff (secondary) |
| --- | ---: | ---: | ---: | ---: |
| UniformNoise(0.0005) / timing-free / ZZZXXX | 0.000496520911 [0.000485021362, 0.000507517852] | 0.993042 | 0.371957 | 0.0005 |
| UniformNoise(0.0005) / timing-free / XXXZZZ | 0.000496941379 [0.00048557676, 0.000507975684] | 0.993883 | 0.372272 | 0.0005 |
| UniformNoise(0.0005) / timing-free / ZXZXZX | 0.000528472418 [0.000517289775, 0.000542212636] | 1.056945 | 0.395892 | 0.000534075843 |
| UniformNoise(0.0005) / timing-free / XZXZXZ | 0.000525387476 [0.000514404965, 0.000538681817] | 1.050775 | 0.393581 | 0.000534075843 |
| UniformNoise(0.001) / timing-free / ZZZXXX | 0.000993711916 [0.00096715895, 0.00101756425] | 0.993712 | 0.744415 | 0.001 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 0.000994388588 [0.000967800198, 0.00101816736] | 0.994389 | 0.744922 | 0.001 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 0.00105757277 [0.00102930567, 0.00110223742] | 1.057573 | 0.792255 | 0.00106815169 |
| UniformNoise(0.001) / timing-free / XZXZXZ | 0.00105437034 [0.00102769623, 0.00109909031] | 1.054370 | 0.789856 | 0.00106815169 |
| UniformNoise(0.0015) / timing-free / ZZZXXX | 0.0015263118 [0.00149388549, 0.0015621277] | 1.017541 | 1.143400 | 0.0015 |
| UniformNoise(0.0015) / timing-free / XXXZZZ | 0.00152681718 [0.00149400936, 0.00156202883] | 1.017878 | 1.143778 | 0.0015 |
| UniformNoise(0.0015) / timing-free / ZXZXZX | 0.00162073867 [0.00159634901, 0.00165244409] | 1.080492 | 1.214137 | 0.00160222753 |
| UniformNoise(0.0015) / timing-free / XZXZXZ | 0.00161110236 [0.0015854094, 0.00164047594] | 1.074068 | 1.206919 | 0.00160222753 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZZZXXX | 0.000360563418 [0.000345337357, 0.000374017211] | 0.721127 | 0.270107 | 0.00035579716 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XXXZZZ | 0.000361184272 [0.000345909003, 0.000374666049] | 0.722369 | 0.270573 | 0.00035579716 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZXZXZX | 0.000377423771 [0.000363355107, 0.000389591416] | 0.754848 | 0.282738 | 0.000371305157 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XZXZXZ | 0.000378336365 [0.000364364236, 0.000390691287] | 0.756673 | 0.283422 | 0.000371305157 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZZZXXX | 0.000710507773 [0.00069299658, 0.000729936658] | 0.710508 | 0.532260 | 0.000711594321 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XXXZZZ | 0.000708913222 [0.000691373233, 0.00072820475] | 0.708913 | 0.531065 | 0.000711594321 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZXZXZX | 0.000749006563 [0.000729571514, 0.000763308439] | 0.749007 | 0.561100 | 0.000742610313 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XZXZXZ | 0.000754886084 [0.000737226081, 0.000768807211] | 0.754886 | 0.565505 | 0.000742610313 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZZZXXX | 0.00143584705 [0.00139151664, 0.00147537813] | 0.717924 | 1.075630 | 0.00142318864 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XXXZZZ | 0.00144366814 [0.00139940188, 0.00148012513] | 0.721834 | 1.081489 | 0.00142318864 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZXZXZX | 0.00150021573 [0.00146129604, 0.0015393221] | 0.750108 | 1.123851 | 0.00148522063 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XZXZXZ | 0.0015047515 [0.00146776411, 0.0015476182] | 0.752376 | 1.127248 | 0.00148522063 |
| H2_LIKE / timing-free / ZZZXXX | 0.00106773642 [0.00103319184, 0.00111453394] | 0.834169 | 0.799869 | 0.00107332623 |
| H2_LIKE / timing-free / XXXZZZ | 0.00104801909 [0.00102414417, 0.00109235873] | 0.818765 | 0.785098 | 0.00107318425 |
| H2_LIKE / timing-free / ZXZXZX | 0.00115892877 [0.00111980014, 0.00118329091] | 0.905413 | 0.868184 | 0.00113618407 |
| H2_LIKE / timing-free / XZXZXZ | 0.0011571009 [0.00111706018, 0.00118220283] | 0.903985 | 0.866814 | 0.0011360424 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 0.00168402644 [0.00164593346, 0.00171519138] | 1.684026 | 1.261548 | 0.0016855869 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 0.00160707193 [0.00158170179, 0.00163476771] | 1.607072 | 1.203899 | 0.00159884911 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 0.00178026283 [0.00174067405, 0.00181353101] | 1.780263 | 1.333641 | 0.00178461589 |
| JONES_MURALI_2026 / SEQ / XZXZXZ | 0.00175760529 [0.00172412473, 0.00179129852] | 1.757605 | 1.316668 | 0.00175750861 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 0.00161963997 [0.00159556139, 0.00165137082] | 1.619640 | 1.213314 | 0.00162061299 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 0.0015722272 [0.0015404283, 0.00159641671] | 1.572227 | 1.177796 | 0.00155981329 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 0.00171858448 [0.00168741803, 0.00175758691] | 1.718584 | 1.287436 | 0.00172898262 |
| JONES_MURALI_2026 / PAR_MR / XZXZXZ | 0.00169304215 [0.00165403923, 0.00172649834] | 1.693042 | 1.268302 | 0.00171026021 |

### Trigger rates and failure share through bare

Rates below are per shot and exclusive by trigger cause; reset + measurement + 2q + idle + mixed = total. Mixed-category fractional and overlapping contributor rates, count errors and raw counts are retained in JSON. Bare failure share is conditional on a logical failure, not on a shot.

| Preset / policy / order | Trigger | Reset | Measurement | 2q | Idle | Mixed | Failure share through bare |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| UniformNoise(0.0005) / timing-free / ZZZXXX | 0.026705396 | 0.005912566 | 0.0059039811 | 0.014874623 | 0 | 1.4226415e-05 | 0.9856058 |
| UniformNoise(0.0005) / timing-free / XXXZZZ | 0.026730151 | 0.005917434 | 0.0059152264 | 0.01488334 | 0 | 1.4150943e-05 | 0.9854335 |
| UniformNoise(0.0005) / timing-free / ZXZXZX | 0.027086679 | 0.0059019057 | 0.0058843585 | 0.015284585 | 0 | 1.5830189e-05 | 0.9877445 |
| UniformNoise(0.0005) / timing-free / XZXZXZ | 0.027081302 | 0.0059057925 | 0.0059130377 | 0.01524766 | 0 | 1.4811321e-05 | 0.9888155 |
| UniformNoise(0.001) / timing-free / ZZZXXX | 0.052745143 | 0.0116925 | 0.011678214 | 0.029315929 | 0 | 5.85e-05 | 0.9847901 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 0.052641214 | 0.011656071 | 0.011629 | 0.029298429 | 0 | 5.7714286e-05 | 0.9872177 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 0.053393929 | 0.011655214 | 0.0116275 | 0.030051929 | 0 | 5.9285714e-05 | 0.9902190 |
| UniformNoise(0.001) / timing-free / XZXZXZ | 0.053325571 | 0.011630143 | 0.011676786 | 0.029957429 | 0 | 6.1214286e-05 | 0.9902668 |
| UniformNoise(0.0015) / timing-free / ZZZXXX | 0.0780675 | 0.017178 | 0.0172745 | 0.043491333 | 0 | 0.00012366667 | 0.9855541 |
| UniformNoise(0.0015) / timing-free / XXXZZZ | 0.078035 | 0.0172175 | 0.0172245 | 0.043455167 | 0 | 0.00013783333 | 0.9865724 |
| UniformNoise(0.0015) / timing-free / ZXZXZX | 0.079196167 | 0.017184833 | 0.017233 | 0.044649 | 0 | 0.00012933333 | 0.9896191 |
| UniformNoise(0.0015) / timing-free / XZXZXZ | 0.079176167 | 0.017298333 | 0.017189167 | 0.044556 | 0 | 0.00013266667 | 0.9884832 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZZZXXX | 0.014992544 | 0 | 0 | 0.014992544 | 0 | 0 | 0.9781752 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XXXZZZ | 0.014977883 | 0 | 0 | 0.014977883 | 0 | 0 | 0.9767418 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZXZXZX | 0.015367097 | 0 | 0 | 0.015367097 | 0 | 0 | 0.9809568 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XZXZXZ | 0.015352932 | 0 | 0 | 0.015352932 | 0 | 0 | 0.9811424 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZZZXXX | 0.029747692 | 0 | 0 | 0.029747692 | 0 | 0 | 0.9743202 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XXXZZZ | 0.029736462 | 0 | 0 | 0.029736462 | 0 | 0 | 0.9755457 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZXZXZX | 0.030435654 | 0 | 0 | 0.030435654 | 0 | 0 | 0.9817895 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XZXZXZ | 0.030462462 | 0 | 0 | 0.030462462 | 0 | 0 | 0.9812802 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZZZXXX | 0.058562714 | 0 | 0 | 0.058562714 | 0 | 0 | 0.9780347 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XXXZZZ | 0.058563714 | 0 | 0 | 0.058563714 | 0 | 0 | 0.9762770 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZXZXZX | 0.059930286 | 0 | 0 | 0.059930286 | 0 | 0 | 0.9773905 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XZXZXZ | 0.059961 | 0 | 0 | 0.059961 | 0 | 0 | 0.9802702 |
| H2_LIKE / timing-free / ZZZXXX | 0.050054833 | 0.00058683333 | 0.011780667 | 0.0376565 | 0 | 3.0833333e-05 | 0.9819724 |
| H2_LIKE / timing-free / XXXZZZ | 0.050067333 | 0.00058283333 | 0.011786083 | 0.03766775 | 0 | 3.0666667e-05 | 0.9825064 |
| H2_LIKE / timing-free / ZXZXZX | 0.050868167 | 0.0005755 | 0.011726667 | 0.0385345 | 0 | 3.15e-05 | 0.9830538 |
| H2_LIKE / timing-free / XZXZXZ | 0.050961333 | 0.00058041667 | 0.0118125 | 0.038535833 | 0 | 3.2583333e-05 | 0.9870682 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 0.104883 | 0.056740833 | 0.011336333 | 0.028508333 | 0.0080425 | 0.000255 | 0.9876224 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 0.1002355 | 0.056944167 | 0.011315167 | 0.028585667 | 0.0031731667 | 0.00021733333 | 0.9703557 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 0.1049325 | 0.057005 | 0.011410833 | 0.029188667 | 0.0070796667 | 0.00024833333 | 0.9865357 |
| JONES_MURALI_2026 / SEQ / XZXZXZ | 0.10309317 | 0.0567525 | 0.011346833 | 0.029264 | 0.0054881667 | 0.00024166667 | 0.9839851 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 0.101841 | 0.056872 | 0.011356833 | 0.028642167 | 0.0047356667 | 0.00023433333 | 0.9901304 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 0.099039 | 0.057050333 | 0.011368 | 0.028569167 | 0.001843 | 0.0002085 | 0.9774322 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 0.102058 | 0.0569235 | 0.011360833 | 0.0293445 | 0.0041976667 | 0.0002315 | 0.9909013 |
| JONES_MURALI_2026 / PAR_MR / XZXZXZ | 0.101037 | 0.056980333 | 0.011336833 | 0.029242167 | 0.0032503333 | 0.00022733333 | 0.9894315 |

### Per-position trigger histogram

All-shot probabilities, not conditional fractions of triggers. Positions refer to the **executed** order; JSON gives each position's actual kind/support label.

| Preset / policy / order | F0 | F1 | F2 | F3 | F4 | F5 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| UniformNoise(0.0005) / timing-free / ZZZXXX | 0.0041272453 | 0.004292 | 0.0043358679 | 0.0049256981 | 0.0045382075 | 0.0044863774 |
| UniformNoise(0.0005) / timing-free / XXXZZZ | 0.0041321887 | 0.004298717 | 0.0043423585 | 0.0049285094 | 0.0045405094 | 0.0044878679 |
| UniformNoise(0.0005) / timing-free / ZXZXZX | 0.0041348491 | 0.0046220943 | 0.0044503019 | 0.0046838302 | 0.0045166226 | 0.0046789811 |
| UniformNoise(0.0005) / timing-free / XZXZXZ | 0.0041142264 | 0.0046154151 | 0.004440434 | 0.0046943019 | 0.0045248679 | 0.0046920566 |
| UniformNoise(0.001) / timing-free / ZZZXXX | 0.0082241429 | 0.0085527143 | 0.0086383571 | 0.0096523571 | 0.0089073571 | 0.0087702143 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 0.0082159286 | 0.0085645714 | 0.0085614286 | 0.0096987143 | 0.0088719286 | 0.0087286429 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 0.0081794286 | 0.0092400714 | 0.0088061429 | 0.0092166429 | 0.0088253571 | 0.0091262857 |
| UniformNoise(0.001) / timing-free / XZXZXZ | 0.0082176429 | 0.0092106429 | 0.0087760714 | 0.0092085 | 0.0088340714 | 0.0090786429 |
| UniformNoise(0.0015) / timing-free / ZZZXXX | 0.0122815 | 0.012761167 | 0.012790667 | 0.014276833 | 0.013123 | 0.012834333 |
| UniformNoise(0.0015) / timing-free / XXXZZZ | 0.0123195 | 0.012802167 | 0.0127465 | 0.014260333 | 0.013085833 | 0.012820667 |
| UniformNoise(0.0015) / timing-free / ZXZXZX | 0.012311 | 0.013661667 | 0.013146 | 0.013634 | 0.012988 | 0.0134555 |
| UniformNoise(0.0015) / timing-free / XZXZXZ | 0.012359333 | 0.013708 | 0.0130715 | 0.013647333 | 0.013053 | 0.013337 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZZZXXX | 0.0021317767 | 0.0023226602 | 0.0023916893 | 0.0029736602 | 0.0026050485 | 0.0025677087 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XXXZZZ | 0.0021303883 | 0.0023251068 | 0.0023816796 | 0.0029757476 | 0.0026043204 | 0.0025606408 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZXZXZX | 0.0021352913 | 0.0026676602 | 0.0024867379 | 0.0027446117 | 0.0025691456 | 0.0027636505 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XZXZXZ | 0.0021318058 | 0.0026521068 | 0.0024873204 | 0.0027527573 | 0.0025647961 | 0.0027641456 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZZZXXX | 0.0042441154 | 0.0046238846 | 0.0047525385 | 0.0058872308 | 0.0051719231 | 0.005068 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XXXZZZ | 0.0042424231 | 0.0046449231 | 0.0047450385 | 0.0058951154 | 0.0051401538 | 0.0050688077 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZXZXZX | 0.0042502308 | 0.0053162692 | 0.0049733462 | 0.0054149231 | 0.0050577692 | 0.0054231154 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XZXZXZ | 0.0042679231 | 0.0053038077 | 0.0049524231 | 0.0054332308 | 0.0050598462 | 0.0054452308 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZZZXXX | 0.0085412857 | 0.0091627143 | 0.0093075714 | 0.011602571 | 0.010093714 | 0.0098548571 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XXXZZZ | 0.0084777143 | 0.0092495714 | 0.0093677143 | 0.011542 | 0.010098286 | 0.0098284286 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZXZXZX | 0.0084355714 | 0.010544857 | 0.0097851429 | 0.010729 | 0.0098934286 | 0.010542286 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XZXZXZ | 0.0084618571 | 0.010473857 | 0.0097755714 | 0.010653286 | 0.0099671429 | 0.010629286 |
| H2_LIKE / timing-free / ZZZXXX | 0.0075174167 | 0.0080315 | 0.0081065833 | 0.0095175 | 0.0085436667 | 0.0083381667 |
| H2_LIKE / timing-free / XXXZZZ | 0.0075200833 | 0.0080280833 | 0.008094 | 0.0095558333 | 0.0085334167 | 0.0083359167 |
| H2_LIKE / timing-free / ZXZXZX | 0.0075953333 | 0.0087826667 | 0.0083366667 | 0.0088705 | 0.0084205 | 0.0088625 |
| H2_LIKE / timing-free / XZXZXZ | 0.0075891667 | 0.0088296667 | 0.0083673333 | 0.00889775 | 0.0084348333 | 0.0088425833 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 0.016258167 | 0.016366667 | 0.0161765 | 0.0200145 | 0.018369 | 0.017698167 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 0.016319333 | 0.017397 | 0.017629 | 0.0170765 | 0.016047 | 0.015766667 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 0.016220667 | 0.018129 | 0.016360333 | 0.019019333 | 0.0160975 | 0.019105667 |
| JONES_MURALI_2026 / SEQ / XZXZXZ | 0.016335 | 0.017001 | 0.018293667 | 0.0166335 | 0.018771 | 0.016059 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 0.016124 | 0.016331167 | 0.016280333 | 0.018769333 | 0.017464167 | 0.016872 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 0.016315833 | 0.0169565 | 0.017065 | 0.0168895 | 0.016033833 | 0.015778333 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 0.016299333 | 0.017649 | 0.016328 | 0.018002833 | 0.015922 | 0.017856833 |
| JONES_MURALI_2026 / PAR_MR / XZXZXZ | 0.016297167 | 0.017011833 | 0.017541167 | 0.016671167 | 0.017562167 | 0.0159535 |

### Logical failure split

Mutually exclusive shares **among failed MC shots**; X-only + Z-only + both XZ = 1. JSON also reports absolute per-shot rates, counts, and inclusive X/Z rates (which overlap through XZ).

| Preset / policy / order | X only | Z only | Both XZ |
| --- | ---: | ---: | ---: |
| UniformNoise(0.0005) / timing-free / ZZZXXX | 0.4187325 | 0.4951020 | 0.0861655 |
| UniformNoise(0.0005) / timing-free / XXXZZZ | 0.4996508 | 0.4126509 | 0.0876983 |
| UniformNoise(0.0005) / timing-free / ZXZXZX | 0.4008692 | 0.4694481 | 0.1296827 |
| UniformNoise(0.0005) / timing-free / XZXZXZ | 0.4744166 | 0.4021136 | 0.1234698 |
| UniformNoise(0.001) / timing-free / ZZZXXX | 0.4123989 | 0.5000963 | 0.0875048 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 0.4975493 | 0.4148006 | 0.0876502 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 0.4116369 | 0.4578666 | 0.1304966 |
| UniformNoise(0.001) / timing-free / XZXZXZ | 0.4709683 | 0.4065279 | 0.1225038 |
| UniformNoise(0.0015) / timing-free / ZZZXXX | 0.4221639 | 0.4865138 | 0.0913224 |
| UniformNoise(0.0015) / timing-free / XXXZZZ | 0.4894498 | 0.4203937 | 0.0901565 |
| UniformNoise(0.0015) / timing-free / ZXZXZX | 0.4061758 | 0.4668778 | 0.1269464 |
| UniformNoise(0.0015) / timing-free / XZXZXZ | 0.4688867 | 0.4019284 | 0.1291849 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZZZXXX | 0.3657674 | 0.5061130 | 0.1281196 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XXXZZZ | 0.5104712 | 0.3701168 | 0.1194120 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / ZXZXZX | 0.3715731 | 0.4715731 | 0.1568537 |
| TWO_QUBIT_ONLY(0.0005) / timing-free / XZXZXZ | 0.4723544 | 0.3744965 | 0.1531490 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZZZXXX | 0.3703927 | 0.5022155 | 0.1273917 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XXXZZZ | 0.5021221 | 0.3729790 | 0.1248989 |
| TWO_QUBIT_ONLY(0.001) / timing-free / ZXZXZX | 0.3772182 | 0.4678993 | 0.1548825 |
| TWO_QUBIT_ONLY(0.001) / timing-free / XZXZXZ | 0.4748425 | 0.3650808 | 0.1600767 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZZZXXX | 0.3685934 | 0.5098266 | 0.1215800 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XXXZZZ | 0.5045915 | 0.3725847 | 0.1228238 |
| TWO_QUBIT_ONLY(0.002) / timing-free / ZXZXZX | 0.3684540 | 0.4803396 | 0.1512064 |
| TWO_QUBIT_ONLY(0.002) / timing-free / XZXZXZ | 0.4657839 | 0.3799325 | 0.1542837 |
| H2_LIKE / timing-free / ZZZXXX | 0.3929432 | 0.4988914 | 0.1081654 |
| H2_LIKE / timing-free / XXXZZZ | 0.5087962 | 0.3915794 | 0.0996244 |
| H2_LIKE / timing-free / ZXZXZX | 0.3898496 | 0.4743213 | 0.1358292 |
| H2_LIKE / timing-free / XZXZXZ | 0.4665857 | 0.3923798 | 0.1410345 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 0.3528841 | 0.5854704 | 0.0616455 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 0.4309199 | 0.5142831 | 0.0547970 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 0.3435953 | 0.5589520 | 0.0974527 |
| JONES_MURALI_2026 / SEQ / XZXZXZ | 0.3915829 | 0.5142644 | 0.0941527 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 0.3701093 | 0.5727882 | 0.0571026 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 0.4394083 | 0.4981035 | 0.0624881 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 0.3676025 | 0.5316121 | 0.1007854 |
| JONES_MURALI_2026 / PAR_MR / XZXZXZ | 0.4140112 | 0.4861489 | 0.0998399 |

Idle exposure: every Jones row contains all seven data qubits' MC-weighted idle seconds and summed Z probabilities, with errors and seven conditional path tables in `results/stageA_check_order.json`. Exposure follows the executed flagged prefix and the permuted bare round. H2/uniform/CX-only rows have no timing assigned.

Assumptions made:
- A1 PAR_MR is retained as an explicit model assumption, using only existing sourced durations. No real hardware schedule, routing, layout, transport events or recooling model is introduced.
- The batched MC samples exactly the original independent per-location Pauli distribution. It has a separate NumPy RNG stream; exhaustive/mixed-fault and tableau oracles validate propagation and decisions. No correction or circuit was changed to obtain any outcome.
- All comparisons use the same preset and timing policy. MC is primary; P2/survival-corrected P2/p_eff_P2 are secondary. Raw P2 determines initial shot allocation only, never the reported p_L or primary p_eff.
- Calibration remains the approved baseline uniform gadget, with lower analysis-only points to bracket all targets. Joint MC uncertainty and threshold ratios use the existing approved bootstrap method. There is no extrapolation and no new pseudothreshold experiment outside one-round uniform MC.
- Significance means the specified per-comparison ≥3σ rule. Null results do not establish equivalence. Descriptive trigger/failure-share changes are not additional significance claims.
- The observed Jones gain supports XXXZZZ for a later Jones comparison under each policy. It does not justify changing the frozen default or claiming the same gain for H2/uniform/CX-only.

Open questions for the human:
- Review the Jones XXXZZZ improvement and the significant losses for alternating orders, including in timing-free noise. The prompt's interleaving/unbiased-noise expectations are not borne out by this unchanged gadget.
- Stop here. Stage B1 (verified ancilla preparation), B2 and B3 have **not started**. No decoder/circuit has been changed, and no Steane-style comparison or paper-matching claim is made.


## Stage A3 — Multi-round memory and alternating-order mechanism

Changes:
- `src/memory_mc.py`: reusable `RoundAdapter`/`RoundStep` driver, arbitrary incoming frames, final endpoint scoring only, current/carry trigger attribution, absorbing-abort interface and count merging. Only the flag adapter is implemented.
- `src/memory_fit.py`: separate marginal binomial fits for epsilon_X/epsilon_Z, observed-curvature MC sampling errors, epsilon_tot with propagated X/Z covariance, goodness of fit, held-out R2 diagnostic and later-window fits.
- `src/memory_experiment.py`: 12 design/preset-policy cases × five R points, cached independent streams, every design-pair contrast and four scientific figures. Union endpoints are descriptive; **no union fit** is included.
- `src/pair_mechanism.py`: exact matched-physical-pair classes, full positive/negative contributions by types/positions/trigger paths, and scalar propagation traces.
- `src/ion_params.py`: user-requested analysis round lengths {1,2,4,8,16}, no hardware values. `src/check_order.py` / `src/order_experiment.py`: old source identity compatibility for validated literal analysis assignments only.
- `tests/test_memory.py`, `tests/fixtures/stageA2_artifacts.json`: carry/tableau/endpoint, single-fault, cache, fit/covariance, R2 diagnostic and pair-conservation checks, plus all 68 approved hashes.
- `results/stageA3_memory.json`, `stageA3_pair_mechanism.json`, `stageA3_validation.json`, `stageA3_cache/`, `figures/memory_*.png`, README and noise-model docs: new A3 checkpoint outputs and methods. No Stage B implementation.

Frozen-result checks:
- All **68 approved Stage 0–A2 result/figure/cache files remain SHA-256 byte-identical**. A1/A2 source identities still match approved values. Existing flag circuits, decoder, timing builders and scalar kernels are unchanged.
- All **119 test nodes** are covered by passing runs: **111 fast tests passed in 58.73 s**; all **8 slow nodes passed in the 115-test full run (414.70 s)**. The focused memory run has **20 passing tests**, including covariance propagation and the R2 trend check. `git diff --check` passes.
- Original flagged FT 0/1218, ignoring-flags 48/1218 and frame/tableau 0/1176 — pass. Uniform malignant count 92840 and all exact A dictionaries — pass; uniform seed 2026 remains exactly 542/3M.
- Jones seed 3326 remains exactly 6121/3M, including trigger metrics and Stage 4 empty replay/pair counts. Approved report reconstruction, idle limits and synthetic transport regressions pass.

FT checks: inject every supported physical payload alone in round 1, 8 or 16 of a 16-round chain and evaluate both endpoint bases. Counts below are **per order**, with all three orders passing each family. Bare-only faults are separate because they do not execute on ideal input without a trigger. Every incoming single data Pauli also traverses 16 noiseless rounds with zero failures in 42 basis trials per design.

| Preset / policy | Flagged failures / basis runs | Bare-only failures / basis runs | Incoming failures / basis runs |
| --- | ---: | ---: | ---: |
| UniformNoise(0.001) / timing-free | 0/3384 | 0/2232 | 0/42 |
| H2_LIKE / timing-free | 0/3456 | 0/2268 | 0/42 |
| JONES_MURALI_2026 / SEQ | 0/3996 | 0/2664 | 0/42 |
| JONES_MURALI_2026 / PAR_MR | 0/3978 | 0/2664 | 0/42 |

Cross-check:
- Six Jones order/policy cases × 90 explicit mixed-fault trajectories × four rounds × both bases = **4,320 tableau round-state/log comparisons** matching scalar and batched frames, syndromes, triggers, corrections and logical decisions.
- Five current/carry trigger components are separately replayed with the scalar primitive through the actual prefix, including cancellations; zero mismatches. R=1 matches frozen A2 on the same injected shots under every preset/policy and each A3 order.
- A toy adapter leaves X0 in round one and adds X1 in round two: final ideal decoding reports failure, while an erroneous intermediate ideal correction would make it pass. Two logical X changes also cancel at R=2, so logical failure is not scored as absorbing first passage. Pre-existing logical frames persist.
- Singletons, incoming errors, generic abort handling, cache/RNG resumption, MLE influence derivatives, covariance propagation and pair-weight conservation pass. No circuit/decoder change was needed.

### Metric definition and sampling

The user's clarification is implemented exactly:

- Fit **epsilon_X** (logical-X marginal / |0_L>) and **epsilon_Z** (logical-Z marginal / |+_L>) separately to p_b(R)=[1-(1-2 epsilon_b)^R]/2 on **R={2,4,8,16}**.
- Primary headline and design-comparison metric: **epsilon_tot=1-(1-epsilon_X)(1-epsilon_Z)**, with propagated uncertainty **including shared-trajectory X/Z covariance**. It is a derived metric, not a fit to union data.
- Union p_L(R) is descriptive at every R, with R=1 retained for continuity. **The union is not fitted.**
- Each marginal has goodness of fit. Compare held-out R=2 with the R≥4 trend; if it is outside its one-sigma MC error bar, report R≥4 too. All R≥4 fits are shown regardless, and R≥8 is an additional sensitivity diagnostic.
- Stage B3 must use this multi-round epsilon_tot as primary, with one-round union p_L secondary. Epsilon_X/Z and fit quality remain required diagnostics.

There are **60 independent R-point samples**, 1M trajectories each: **60M endpoints / 372M noisy rounds**. The full output Pauli frame becomes the next round's input; there is no ideal cleanup or added inter-round time. Actual noisy protocol corrections occur each round; one ideal decoder scores the endpoint. H2 uses timing-free Stage 1 channels; no timing is invented. Jones runs all three orders under both SEQ and PAR_MR, nominal parameters and empty transport.

Sampling stderr for each marginal uses its observed likelihood curvature and endpoint-binomial variance. Within the same R, Cov(k_X,k_Z)=N(p_XZ-p_X p_Z), propagated through both MLE influence functions. Then Var(epsilon_tot) includes the cross term 2(1-epsilon_X)(1-epsilon_Z) Cov(epsilon_X,epsilon_Z). Across designs and R points the streams are independent. Expected-Fisher errors are secondary checks, and Pearson inflation is only diagnostic. Maximum relative sampling stderr across primary marginal/combined metrics is below **0.90%**, satisfying ≤3%. Positive fitted X/Z estimator correlation is 0.108–0.245; retaining covariance increases headline stderr by **5.2–11.6%** over the independence approximation.

### Prediction test: epsilon_Z and epsilon_tot

**Near-equality of grouped orders is not established over R≤16.** Under Jones, XXXZZZ retains lower epsilon_Z by 12.31% (SEQ) / 11.68% (PAR_MR), and lower epsilon_tot by 4.55% / 3.50%. The totals are closer than their one-round product-summary counterparts, but epsilon_Z does not converge near zero difference in this fit window. Later-window results are explicitly reported below. The grouped check lists are cyclic shifts when all checks execute; the actual adaptive protocol can terminate a flagged prefix, run bare checks and restart flagged checks at a boundary, so cyclic equivalence of the full stochastic process is not assumed.

Fit quality limits stationary interpretation: **13 of 24 primary marginal fits have GOF p<0.05**, including seven with p<0.0027. Some R≥4 fits also remain poor. These are model diagnostics, not code/FT failures. Deltas and sigmas quantify finite-window fitted metrics; the data do not establish an infinite-R limit.

For a consistent boundary comparison, the R1 product summary below uses 1-(1-p_X(1))(1-p_Z(1)) from approved A2 marginal counts; it is not a union fit. Prior one-round union differences are also retained for continuity.

| Jones policy | Metric | R1 corresponding marginal/product relative difference | Memory relative difference [95% sampling CI] | Δepsilon ± combined stderr | σ | R≥4 relative difference |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| SEQ | epsilon_Z | -20.8026% | -12.3134% [-13.5000%, -11.1267%] | -0.0001895288183 ± 9.953e-06 | -19.0424 | -11.8038% |
| SEQ | epsilon_tot | -10.5204% | -4.5524% [-5.5963%, -3.5086%] | -0.0001130520662 ± 1.35377e-05 | -8.3509 | -4.1967% |
| PAR_MR | epsilon_Z | -17.2916% | -11.6820% [-12.9489%, -10.4152%] | -0.0001589444678 ± 9.3587e-06 | -16.9836 | -11.6580% |
| PAR_MR | epsilon_tot | -6.5924% | -3.4970% [-4.5935%, -2.4004%] | -8.022579968e-05 ± 1.30659e-05 | -6.1401 | -3.5815% |

The earlier A2 union differences were -9.9426% SEQ and -7.0673% PAR_MR; those are historical one-round descriptive numbers, not the definition of the new combined metric. Uniform/H2 epsilon_tot grouped-order comparisons remain below 3σ; ZXZXZX raises epsilon_tot under every family. Epsilon_X moves in the opposite direction to epsilon_Z when the grouped order is reversed, so a small combined difference alone does not imply marginal convergence.

### Primary epsilon_tot, R1 continuity, and R≥4 result

All errors below are propagated sampling stderr. Epsilon_tot has no standalone decay-curve goodness-of-fit statistic: the required quality belongs to its two marginal fits, reported in the next table. R1 is descriptive union p_L.

| Preset / policy / order | Union p_L(1) ± stderr | epsilon_tot ± stderr | Relative stderr | X/Z estimator correlation | R≥4 epsilon_tot ± stderr |
| --- | ---: | ---: | ---: | ---: | ---: |
| UniformNoise(0.001) / timing-free / ZZZXXX | 0.000739 ± 2.7175e-05 | 0.0008951398211 ± 5.9139e-06 | 0.6607% | 0.16286 | 0.0008959284581 ± 6.1259e-06 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 0.000735 ± 2.7101e-05 | 0.0009094372086 ± 5.9571e-06 | 0.6550% | 0.16136 | 0.0009134215839 ± 6.1831e-06 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 0.000889 ± 2.9803e-05 | 0.001040955505 ± 6.5466e-06 | 0.6289% | 0.22389 | 0.00104266555 ± 6.7864e-06 |
| H2_LIKE / timing-free / ZZZXXX | 0.000836 ± 2.8902e-05 | 0.001055299128 ± 6.5081e-06 | 0.6167% | 0.19333 | 0.001057177067 ± 6.7461e-06 |
| H2_LIKE / timing-free / XXXZZZ | 0.000794 ± 2.8167e-05 | 0.001053202721 ± 6.4851e-06 | 0.6157% | 0.18717 | 0.001058364265 ± 6.7305e-06 |
| H2_LIKE / timing-free / ZXZXZX | 0.000927 ± 3.0433e-05 | 0.001199683679 ± 7.0928e-06 | 0.5912% | 0.24500 | 0.001200069269 ± 7.3457e-06 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 0.002037 ± 4.5087e-05 | 0.002483326068 ± 9.6893e-06 | 0.3902% | 0.11061 | 0.002489353591 ± 1.0051e-05 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 0.001819 ± 4.2611e-05 | 0.002370274002 ± 9.4544e-06 | 0.3989% | 0.10809 | 0.002384883362 ± 9.8234e-06 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 0.002243 ± 4.7307e-05 | 0.002836071493 ± 1.0654e-05 | 0.3756% | 0.17243 | 0.002848657861 ± 1.1068e-05 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 0.001883 ± 4.3353e-05 | 0.00229414475 ± 9.326e-06 | 0.4065% | 0.11510 | 0.002306870085 ± 9.6898e-06 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 0.001768 ± 4.201e-05 | 0.00221391895 ± 9.1511e-06 | 0.4133% | 0.11257 | 0.002224250022 ± 9.501e-06 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 0.002145 ± 4.6264e-05 | 0.002633287736 ± 1.0268e-05 | 0.3899% | 0.17415 | 0.002635485343 ± 1.0633e-05 |

### Marginal fits and goodness of fit (R≥2)

Chi-square has three degrees of freedom; the p-value and binomial deviance are also stored with standardized per-point residuals in JSON. Marginal fits are separate likelihood fits; their estimators are not assumed statistically independent.

| Preset / policy / order | epsilon_X ± stderr | X chi²/3, GOF p | epsilon_Z ± stderr | Z chi²/3, GOF p |
| --- | ---: | ---: | ---: | ---: |
| UniformNoise(0.001) / timing-free / ZZZXXX | 0.0004228076906 ± 3.7698e-06 | 0.41312/3, 0.93752 | 0.0004725319206 ± 3.9873e-06 | 0.17539/3, 0.98146 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 0.0004838091815 ± 4.0347e-06 | 14.156/3, 0.0026999 | 0.0004258340495 ± 3.7833e-06 | 1.0867/3, 0.78029 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 0.0004938097111 ± 4.0768e-06 | 3.2998/3, 0.34768 | 0.0005474161132 ± 4.2945e-06 | 7.0204/3, 0.071251 |
| H2_LIKE / timing-free / ZZZXXX | 0.0004817218839 ± 4.026e-06 | 5.3248/3, 0.1495 | 0.0005738536824 ± 4.3981e-06 | 7.7983/3, 0.050369 |
| H2_LIKE / timing-free / XXXZZZ | 0.0005684571717 ± 4.3772e-06 | 13.205/3, 0.0042132 | 0.0004850212626 ± 4.0399e-06 | 9.9704/3, 0.018819 |
| H2_LIKE / timing-free / ZXZXZX | 0.0005572613086 ± 4.3336e-06 | 4.022/3, 0.2591 | 0.0006427805673 ± 4.658e-06 | 4.6842/3, 0.19644 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 0.0009455690772 ± 5.6665e-06 | 9.5511/3, 0.022794 | 0.001539212423 ± 7.2718e-06 | 12.522/3, 0.0057928 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 0.001021969733 ± 5.895e-06 | 18.919/3, 0.00028411 | 0.001349683604 ± 6.7958e-06 | 50.005/3, 7.9689e-11 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 0.001136013493 ± 6.2224e-06 | 13.231/3, 0.0041618 | 0.001701991486 ± 7.6589e-06 | 16.967/3, 0.0007179 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 0.0009348293304 ± 5.6336e-06 | 15.561/3, 0.0013951 | 0.001360587337 ± 6.8249e-06 | 22.796/3, 4.4541e-05 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 0.001013493939 ± 5.8701e-06 | 13.068/3, 0.0044908 | 0.001201642869 ± 6.4036e-06 | 21.333/3, 8.9755e-05 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 0.001107460887 ± 6.1422e-06 | 4.6595/3, 0.1985 | 0.001527518516 ± 7.244e-06 | 3.8031/3, 0.28353 |

### Held-out R=2 diagnostic and R≥4 marginal fits

Compare R=2 with the independent R≥4 trend. “Outside” means beyond the R2 point's reported **one-sigma MC error bar**; combined sigma also propagates trend-fit uncertainty. Twelve of 24 marginals meet this criterion, covering eight of 12 designs. Every required R≥4 fit is reported, and the remaining ones are shown too. The R≥4 fit has two degrees of freedom. R≥8 sensitivities are retained in JSON without replacing the requested primary R≥2 fit.

| Preset / policy / order | Marginal | R2 deviation / own error bar | Combined trend σ | Outside R2 bar? | R≥4 epsilon ± stderr | R≥4 chi²/2, GOF p |
| --- | --- | ---: | ---: | --- | ---: | ---: |
| UniformNoise(0.001) / timing-free / ZZZXXX | epsilon_X | -0.4181 | -0.4036 | no | 0.000423213522 ± 3.9051e-06 | 0.252/2, 0.88162 |
| UniformNoise(0.001) / timing-free / ZZZXXX | epsilon_Z | -0.3726 | -0.3597 | no | 0.0004729150801 ± 4.1302e-06 | 0.047293/2, 0.97663 |
| UniformNoise(0.001) / timing-free / XXXZZZ | epsilon_X | -3.1132 | -2.9961 | yes | 0.0004869190645 ± 4.1913e-06 | 5.8835/2, 0.052772 |
| UniformNoise(0.001) / timing-free / XXXZZZ | epsilon_Z | -0.9064 | -0.8745 | no | 0.0004267102928 ± 3.9213e-06 | 0.34143/2, 0.84306 |
| UniformNoise(0.001) / timing-free / ZXZXZX | epsilon_X | -0.9744 | -0.9401 | no | 0.0004948254219 ± 4.2257e-06 | 2.4344/2, 0.29606 |
| UniformNoise(0.001) / timing-free / ZXZXZX | epsilon_Z | -0.6293 | -0.6073 | no | 0.0005481113471 ± 4.4498e-06 | 6.6492/2, 0.035987 |
| H2_LIKE / timing-free / ZZZXXX | epsilon_X | -1.3755 | -1.3264 | yes | 0.0004831288952 ± 4.1749e-06 | 3.621/2, 0.16357 |
| H2_LIKE / timing-free / ZZZXXX | epsilon_Z | -0.4157 | -0.4013 | no | 0.000574325645 ± 4.5562e-06 | 7.6327/2, 0.022007 |
| H2_LIKE / timing-free / XXXZZZ | epsilon_X | -3.8700 | -3.7222 | yes | 0.0005726259993 ± 4.5496e-06 | 0.64815/2, 0.72319 |
| H2_LIKE / timing-free / XXXZZZ | epsilon_Z | -0.9636 | -0.9296 | no | 0.0004860165716 ± 4.1874e-06 | 9.1104/2, 0.010512 |
| H2_LIKE / timing-free / ZXZXZX | epsilon_X | 0.2275 | 0.2197 | no | 0.0005570046589 ± 4.4863e-06 | 3.9753/2, 0.13702 |
| H2_LIKE / timing-free / ZXZXZX | epsilon_Z | -0.5346 | -0.5160 | no | 0.0006434229999 ± 4.826e-06 | 4.417/2, 0.10987 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | epsilon_X | -3.3106 | -3.1878 | yes | 0.0009502864242 ± 5.8843e-06 | 0.03267/2, 0.9838 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | epsilon_Z | -0.6962 | -0.6716 | no | 0.001540531113 ± 7.5383e-06 | 12.066/2, 0.0023986 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | epsilon_X | -3.8025 | -3.6602 | yes | 0.001027590165 ± 6.1235e-06 | 6.4217/2, 0.040323 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | epsilon_Z | -5.3125 | -5.1091 | yes | 0.001358689373 ± 7.0655e-06 | 25.924/2, 2.3476e-06 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | epsilon_X | -3.3660 | -3.2414 | yes | 0.001141303703 ± 6.4614e-06 | 3.3284/2, 0.18934 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | epsilon_Z | -3.7490 | -3.6095 | yes | 0.001709304994 ± 7.955e-06 | 4.6176/2, 0.099378 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | epsilon_X | -3.5065 | -3.3758 | yes | 0.0009397847839 ± 5.851e-06 | 4.9064/2, 0.086018 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | epsilon_Z | -4.5402 | -4.3687 | yes | 0.001368371276 ± 7.0923e-06 | 5.0544/2, 0.079881 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | epsilon_X | -2.0958 | -2.0201 | yes | 0.001016632566 ± 6.0901e-06 | 9.12/2, 0.010462 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | epsilon_Z | -4.4945 | -4.3245 | yes | 0.001208846408 ± 6.6545e-06 | 4.0235/2, 0.13375 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | epsilon_X | -0.6775 | -0.6537 | no | 0.001108539425 ± 6.3657e-06 | 4.2333/2, 0.12043 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | epsilon_Z | -0.5943 | -0.5733 | no | 0.001528640477 ± 7.509e-06 | 3.4749/2, 0.17597 |

### Every design pair: primary epsilon_tot and prediction metric epsilon_Z

Direction is second minus first under the same preset/policy; combined uncertainty uses independent design streams, after the X/Z covariance has already been propagated inside each epsilon_tot. Three-sigma labels concern sampling comparisons of finite-window metrics, not proof of stationarity. No family-wise multiplicity correction is claimed.

| Family | First → second | Δepsilon_tot ± stderr | σ_tot | R≥4 Δepsilon_tot ± stderr (σ) | Δepsilon_Z ± stderr | σ_Z | R≥4 Δepsilon_Z ± stderr (σ) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| UniformNoise(0.001) / timing-free | ZZZXXX → XXXZZZ | 1.42973875e-05 ± 8.3941e-06 | 1.7033 | 1.74931259e-05 ± 8.7039e-06 (2.010) | -4.66978711e-05 ± 5.4965e-06 | -8.4959 | -4.62047873e-05 ± 5.6952e-06 (-8.113) |
| UniformNoise(0.001) / timing-free | ZZZXXX → ZXZXZX | 0.000145815684 ± 8.8223e-06 | 16.5281 | 0.000146737091 ± 9.1423e-06 (16.050) | 7.48841926e-05 ± 5.8601e-06 | 12.7786 | 7.5196267e-05 ± 6.0712e-06 (12.386) |
| UniformNoise(0.001) / timing-free | XXXZZZ → ZXZXZX | 0.000131518296 ± 8.8513e-06 | 14.8587 | 0.000129243966 ± 9.1807e-06 (14.078) | 0.000121582064 ± 5.7233e-06 | 21.2432 | 0.000121401054 ± 5.931e-06 (20.469) |
| H2_LIKE / timing-free | ZZZXXX → XXXZZZ | -2.09640782e-06 ± 9.1875e-06 | -0.2282 | 1.18719832e-06 ± 9.5294e-06 (0.125) | -8.88324198e-05 ± 5.9719e-06 | -14.8750 | -8.83090734e-05 ± 6.1881e-06 (-14.271) |
| H2_LIKE / timing-free | ZZZXXX → ZXZXZX | 0.000144384551 ± 9.6262e-06 | 14.9992 | 0.000142892202 ± 9.9735e-06 (14.327) | 6.8926885e-05 ± 6.4063e-06 | 10.7593 | 6.90973549e-05 ± 6.6369e-06 (10.411) |
| H2_LIKE / timing-free | XXXZZZ → ZXZXZX | 0.000146480959 ± 9.6106e-06 | 15.2416 | 0.000141705004 ± 9.9629e-06 (14.223) | 0.000157759305 ± 6.1658e-06 | 25.5861 | 0.000157406428 ± 6.3894e-06 (24.636) |
| JONES_MURALI_2026 / SEQ | ZZZXXX → XXXZZZ | -0.000113052066 ± 1.3538e-05 | -8.3509 | -0.000104470229 ± 1.4054e-05 (-7.434) | -0.000189528818 ± 9.953e-06 | -19.0424 | -0.00018184174 ± 1.0332e-05 (-17.600) |
| JONES_MURALI_2026 / SEQ | ZZZXXX → ZXZXZX | 0.000352745425 ± 1.4401e-05 | 24.4949 | 0.00035930427 ± 1.4951e-05 (24.033) | 0.000162779063 ± 1.0561e-05 | 15.4131 | 0.000168773881 ± 1.0959e-05 (15.400) |
| JONES_MURALI_2026 / SEQ | XXXZZZ → ZXZXZX | 0.000465797491 ± 1.4244e-05 | 32.7018 | 0.000463774499 ± 1.4799e-05 (31.339) | 0.000352307882 ± 1.0239e-05 | 34.4077 | 0.000350615621 ± 1.064e-05 (32.954) |
| JONES_MURALI_2026 / PAR_MR | ZZZXXX → XXXZZZ | -8.02257997e-05 ± 1.3066e-05 | -6.1401 | -8.26200633e-05 ± 1.3571e-05 (-6.088) | -0.000158944468 ± 9.3587e-06 | -16.9836 | -0.000159524868 ± 9.7254e-06 (-16.403) |
| JONES_MURALI_2026 / PAR_MR | ZZZXXX → ZXZXZX | 0.000339142986 ± 1.3871e-05 | 24.4500 | 0.000328615258 ± 1.4386e-05 (22.843) | 0.00016693118 ± 9.9526e-06 | 16.7726 | 0.000160269201 ± 1.0329e-05 (15.517) |
| JONES_MURALI_2026 / PAR_MR | XXXZZZ → ZXZXZX | 0.000419368786 ± 1.3754e-05 | 30.4909 | 0.000411235321 ± 1.4259e-05 (28.840) | 0.000325875647 ± 9.6686e-06 | 33.7045 | 0.000319794068 ± 1.0033e-05 (31.873) |

| Family | First → second | Δepsilon_X ± combined stderr | σ_X |
| --- | --- | ---: | ---: |
| UniformNoise(0.001) / timing-free | ZZZXXX → XXXZZZ | 6.10014909e-05 ± 5.5218e-06 | 11.0474 |
| UniformNoise(0.001) / timing-free | ZZZXXX → ZXZXZX | 7.10020205e-05 ± 5.5526e-06 | 12.7872 |
| UniformNoise(0.001) / timing-free | XXXZZZ → ZXZXZX | 1.00005296e-05 ± 5.7358e-06 | 1.7435 |
| H2_LIKE / timing-free | ZZZXXX → XXXZZZ | 8.67352879e-05 ± 5.9472e-06 | 14.5843 |
| H2_LIKE / timing-free | ZZZXXX → ZXZXZX | 7.55394247e-05 ± 5.9152e-06 | 12.7705 |
| H2_LIKE / timing-free | XXXZZZ → ZXZXZX | -1.11958632e-05 ± 6.1596e-06 | -1.8176 |
| JONES_MURALI_2026 / SEQ | ZZZXXX → XXXZZZ | 7.64006562e-05 ± 8.1768e-06 | 9.3436 |
| JONES_MURALI_2026 / SEQ | ZZZXXX → ZXZXZX | 0.000190444416 ± 8.4159e-06 | 22.6291 |
| JONES_MURALI_2026 / SEQ | XXXZZZ → ZXZXZX | 0.000114043759 ± 8.5714e-06 | 13.3051 |
| JONES_MURALI_2026 / PAR_MR | ZZZXXX → XXXZZZ | 7.86646089e-05 ± 8.1361e-06 | 9.6686 |
| JONES_MURALI_2026 / PAR_MR | ZZZXXX → ZXZXZX | 0.000172631557 ± 8.3345e-06 | 20.7130 |
| JONES_MURALI_2026 / PAR_MR | XXXZZZ → ZXZXZX | 9.39669477e-05 ± 8.4962e-06 | 11.0599 |

### Endpoint Monte Carlo, all R: union descriptive only

Each row has 1M independent trajectories. There is no intermediate ideal correction, no union fit and no absorbing first-logical-failure count. X/Z marginals can fail together; their overlap is retained in the covariance calculation. The generic abort hook is unused in this flag study.

| Preset / policy / order | R | Union p_L ± stderr | X / 0_L marginal ± stderr | Z / +_L marginal ± stderr |
| --- | ---: | ---: | ---: | ---: |
| UniformNoise(0.001) / timing-free / ZZZXXX | 1 | 0.000739 ± 2.7175e-05 | 0.000364 ± 1.9075e-05 | 0.000425 ± 2.0611e-05 |
| UniformNoise(0.001) / timing-free / ZZZXXX | 2 | 0.001623 ± 4.0254e-05 | 0.000834 ± 2.8867e-05 | 0.000934 ± 3.0547e-05 |
| UniformNoise(0.001) / timing-free / ZZZXXX | 4 | 0.00329 ± 5.7264e-05 | 0.001683 ± 4.099e-05 | 0.001896 ± 4.3502e-05 |
| UniformNoise(0.001) / timing-free / ZZZXXX | 8 | 0.006569 ± 8.0783e-05 | 0.0034 ± 5.821e-05 | 0.003762 ± 6.122e-05 |
| UniformNoise(0.001) / timing-free / ZZZXXX | 16 | 0.013024 ± 0.00011338 | 0.006712 ± 8.1651e-05 | 0.007515 ± 8.6363e-05 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 1 | 0.000735 ± 2.7101e-05 | 0.000413 ± 2.0318e-05 | 0.000388 ± 1.9694e-05 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 2 | 0.001576 ± 3.9668e-05 | 0.000881 ± 2.9669e-05 | 0.000827 ± 2.8746e-05 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 4 | 0.003286 ± 5.7229e-05 | 0.001906 ± 4.3616e-05 | 0.001697 ± 4.116e-05 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 8 | 0.006585 ± 8.088e-05 | 0.003782 ± 6.1382e-05 | 0.00338 ± 5.8039e-05 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 16 | 0.013468 ± 0.00011527 | 0.007873 ± 8.838e-05 | 0.006815 ± 8.2271e-05 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 1 | 0.000889 ± 2.9803e-05 | 0.000463 ± 2.1512e-05 | 0.000536 ± 2.3145e-05 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 2 | 0.001816 ± 4.2576e-05 | 0.000959 ± 3.0953e-05 | 0.001075 ± 3.277e-05 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 4 | 0.00361 ± 5.9975e-05 | 0.001943 ± 4.4037e-05 | 0.002142 ± 4.6232e-05 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 8 | 0.007232 ± 8.4733e-05 | 0.003888 ± 6.2232e-05 | 0.004258 ± 6.5114e-05 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 16 | 0.014866 ± 0.00012102 | 0.007949 ± 8.8802e-05 | 0.008855 ± 9.3683e-05 |
| H2_LIKE / timing-free / ZZZXXX | 1 | 0.000836 ± 2.8902e-05 | 0.000407 ± 2.017e-05 | 0.000504 ± 2.2444e-05 |
| H2_LIKE / timing-free / ZZZXXX | 2 | 0.001865 ± 4.3145e-05 | 0.000924 ± 3.0383e-05 | 0.001134 ± 3.3656e-05 |
| H2_LIKE / timing-free / ZZZXXX | 4 | 0.003689 ± 6.0625e-05 | 0.00187 ± 4.3203e-05 | 0.002193 ± 4.6778e-05 |
| H2_LIKE / timing-free / ZZZXXX | 8 | 0.007533 ± 8.6465e-05 | 0.003811 ± 6.1616e-05 | 0.004523 ± 6.7101e-05 |
| H2_LIKE / timing-free / ZZZXXX | 16 | 0.015298 ± 0.00012274 | 0.007775 ± 8.7833e-05 | 0.009264 ± 9.5803e-05 |
| H2_LIKE / timing-free / XXXZZZ | 1 | 0.000794 ± 2.8167e-05 | 0.000478 ± 2.1858e-05 | 0.000395 ± 1.9871e-05 |
| H2_LIKE / timing-free / XXXZZZ | 2 | 0.001772 ± 4.2058e-05 | 0.001021 ± 3.1937e-05 | 0.000942 ± 3.0678e-05 |
| H2_LIKE / timing-free / XXXZZZ | 4 | 0.003704 ± 6.0748e-05 | 0.002251 ± 4.7391e-05 | 0.001841 ± 4.2867e-05 |
| H2_LIKE / timing-free / XXXZZZ | 8 | 0.007612 ± 8.6914e-05 | 0.004575 ± 6.7484e-05 | 0.00382 ± 6.1688e-05 |
| H2_LIKE / timing-free / XXXZZZ | 16 | 0.015335 ± 0.00012288 | 0.009107 ± 9.4995e-05 | 0.007875 ± 8.8391e-05 |
| H2_LIKE / timing-free / ZXZXZX | 1 | 0.000927 ± 3.0433e-05 | 0.0005 ± 2.2355e-05 | 0.000542 ± 2.3275e-05 |
| H2_LIKE / timing-free / ZXZXZX | 2 | 0.002094 ± 4.5712e-05 | 0.001121 ± 3.3463e-05 | 0.001267 ± 3.5572e-05 |
| H2_LIKE / timing-free / ZXZXZX | 4 | 0.00406 ± 6.3589e-05 | 0.002148 ± 4.6297e-05 | 0.002495 ± 4.9888e-05 |
| H2_LIKE / timing-free / ZXZXZX | 8 | 0.008383 ± 9.1174e-05 | 0.004516 ± 6.7049e-05 | 0.005068 ± 7.1009e-05 |
| H2_LIKE / timing-free / ZXZXZX | 16 | 0.016778 ± 0.00012844 | 0.008837 ± 9.3589e-05 | 0.010326 ± 0.00010109 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 1 | 0.002037 ± 4.5087e-05 | 0.000826 ± 2.8728e-05 | 0.001336 ± 3.6527e-05 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 2 | 0.004533 ± 6.7175e-05 | 0.00176 ± 4.1915e-05 | 0.003038 ± 5.5034e-05 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 4 | 0.009125 ± 9.5088e-05 | 0.003784 ± 6.1398e-05 | 0.005907 ± 7.663e-05 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 8 | 0.018562 ± 0.00013497 | 0.007544 ± 8.6528e-05 | 0.012137 ± 0.0001095 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 16 | 0.036941 ± 0.00018862 | 0.015004 ± 0.00012157 | 0.024369 ± 0.00015419 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 1 | 0.001819 ± 4.2611e-05 | 0.000873 ± 2.9534e-05 | 0.001051 ± 3.2402e-05 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 2 | 0.004085 ± 6.3783e-05 | 0.001888 ± 4.341e-05 | 0.002451 ± 4.9447e-05 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 4 | 0.00858 ± 9.223e-05 | 0.004003 ± 6.3143e-05 | 0.005094 ± 7.119e-05 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 8 | 0.017708 ± 0.00013189 | 0.00805 ± 8.936e-05 | 0.01071 ± 0.00010293 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 16 | 0.035719 ± 0.00018559 | 0.016397 ± 0.000127 | 0.021677 ± 0.00014563 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 1 | 0.002243 ± 4.7307e-05 | 0.00102 ± 3.1921e-05 | 0.001441 ± 3.7933e-05 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 2 | 0.004902 ± 6.9842e-05 | 0.002125 ± 4.6049e-05 | 0.003201 ± 5.6487e-05 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 4 | 0.010211 ± 0.00010053 | 0.004513 ± 6.7027e-05 | 0.006696 ± 8.1555e-05 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 8 | 0.020376 ± 0.00014128 | 0.008937 ± 9.4112e-05 | 0.013395 ± 0.00011496 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 16 | 0.04075 ± 0.00019771 | 0.018109 ± 0.00013335 | 0.026882 ± 0.00016174 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 1 | 0.001883 ± 4.3353e-05 | 0.000776 ± 2.7846e-05 | 0.001244 ± 3.5248e-05 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 2 | 0.004 ± 6.3119e-05 | 0.001732 ± 4.1581e-05 | 0.002506 ± 4.9997e-05 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 4 | 0.008393 ± 9.1228e-05 | 0.003637 ± 6.0198e-05 | 0.005309 ± 7.2669e-05 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 8 | 0.0174 ± 0.00013076 | 0.007577 ± 8.6716e-05 | 0.010961 ± 0.00010412 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 16 | 0.034009 ± 0.00018125 | 0.01483 ± 0.00012087 | 0.021474 ± 0.00014496 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 1 | 0.001768 ± 4.201e-05 | 0.000876 ± 2.9584e-05 | 0.000983 ± 3.1337e-05 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 2 | 0.003896 ± 6.2296e-05 | 0.001939 ± 4.3991e-05 | 0.002204 ± 4.6895e-05 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 4 | 0.008112 ± 8.9701e-05 | 0.003885 ± 6.2209e-05 | 0.004713 ± 6.8489e-05 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 8 | 0.016606 ± 0.00012779 | 0.008066 ± 8.9448e-05 | 0.009533 ± 9.7171e-05 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 16 | 0.033033 ± 0.00017872 | 0.016199 ± 0.00012624 | 0.019156 ± 0.00013707 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 1 | 0.002145 ± 4.6264e-05 | 0.001004 ± 3.167e-05 | 0.001347 ± 3.6677e-05 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 2 | 0.004675 ± 6.8214e-05 | 0.002183 ± 4.6672e-05 | 0.00302 ± 5.4871e-05 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 4 | 0.00938 ± 9.6395e-05 | 0.004318 ± 6.5569e-05 | 0.005971 ± 7.7041e-05 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 8 | 0.019122 ± 0.00013695 | 0.008741 ± 9.3084e-05 | 0.012222 ± 0.00010988 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 16 | 0.037616 ± 0.00019027 | 0.017605 ± 0.00013151 | 0.023898 ± 0.00015273 |

### Trigger and bare-round metrics

R=16 summary below; JSON provides every R, every round index, positional histograms, exclusive causes, overlapping contributors, fractional credits and raw counts. Rates are per executed round; `carry` is incoming residual-frame signal, not a current location. Any-bare probability is per whole trajectory. Memory failure share through any bare round and through the last bare round are separate.

| Preset / policy / order | Trigger/round | Reset | Measurement | 2q | Idle | Carry | Mixed | Any-bare probability | Failure share any bare | Failure share last bare |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| UniformNoise(0.001) / timing-free / ZZZXXX | 0.05772944 | 0.01161937 | 0.01157725 | 0.02920556 | 0 | 0.005255125 | 7.2125e-05 | 0.605539 | 0.9996929 | 0.1628532 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 0.05775913 | 0.01162206 | 0.01157394 | 0.02922687 | 0 | 0.0052615 | 7.475e-05 | 0.605055 | 0.9998515 | 0.1647609 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 0.05795031 | 0.01155187 | 0.01164425 | 0.03003531 | 0 | 0.004646687 | 7.21875e-05 | 0.606022 | 0.9996637 | 0.1652092 |
| H2_LIKE / timing-free / ZZZXXX | 0.05597837 | 0.000581625 | 0.01171594 | 0.03741738 | 0 | 0.006216562 | 4.6875e-05 | 0.594033 | 0.9996732 | 0.1615244 |
| H2_LIKE / timing-free / XXXZZZ | 0.05605456 | 0.0005729375 | 0.01176231 | 0.0374495 | 0 | 0.006221875 | 4.79375e-05 | 0.59441 | 0.9994131 | 0.1627649 |
| H2_LIKE / timing-free / ZXZXZX | 0.05616213 | 0.000578 | 0.01168219 | 0.03841825 | 0 | 0.00544025 | 4.34375e-05 | 0.594617 | 0.9994636 | 0.1630707 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 0.1153734 | 0.05639437 | 0.01126012 | 0.02835987 | 0.007924125 | 0.01111112 | 0.00032375 | 0.846203 | 0.9997834 | 0.2127446 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 0.1148024 | 0.05613488 | 0.01120888 | 0.02817569 | 0.00312175 | 0.01583794 | 0.0003233125 | 0.846234 | 0.9996081 | 0.2138638 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 0.1154029 | 0.05635431 | 0.01123769 | 0.02899556 | 0.006944625 | 0.01154481 | 0.0003259375 | 0.846996 | 0.9997546 | 0.2123190 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 0.1111816 | 0.05649688 | 0.01125531 | 0.02834288 | 0.00471775 | 0.01007525 | 0.0002935 | 0.83523 | 0.9999412 | 0.2174130 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 0.1111188 | 0.05636581 | 0.01125856 | 0.02828931 | 0.001812438 | 0.01309531 | 0.0002973125 | 0.835113 | 0.9997881 | 0.2149669 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 0.1115096 | 0.05655113 | 0.01129156 | 0.02901431 | 0.004143625 | 0.01021681 | 0.000292125 | 0.835768 | 0.9998937 | 0.2124894 |

All Jones rows include actual-path mean idle seconds and summed Z probability for every data qubit, normalized per executed round. Independent-round exposure stderr is omitted because round paths within one trajectory can correlate. H2 remains timing-free.

### One-round P2 and survival correction (secondary only)

These are the unchanged one-round A2 pair tables, compared to fresh A3 R=1 MC. They are **not** multi-round P2 estimates or epsilon estimates. Raw and corrected values and both ratios are reported for every design/preset. No P2×R approximation is substituted for the memory simulation.

| Preset / policy / order | Raw P2 | Exact survival-corrected P2 | Raw/R1 MC | Corrected/R1 MC |
| --- | ---: | ---: | ---: | ---: |
| UniformNoise(0.001) / timing-free / ZZZXXX | 0.00076791111111 | 0.00069898337247 | 1.039122 | 0.945850 |
| UniformNoise(0.001) / timing-free / XXXZZZ | 0.00076791111111 | 0.00069898337247 | 1.044777 | 0.950998 |
| UniformNoise(0.001) / timing-free / ZXZXZX | 0.00087614666666 | 0.0007975036993 | 0.985542 | 0.897080 |
| H2_LIKE / timing-free / ZZZXXX | 0.00088465601772 | 0.00080565603442 | 1.058201 | 0.963703 |
| H2_LIKE / timing-free / XXXZZZ | 0.00088442199291 | 0.00080544293281 | 1.113882 | 1.014412 |
| H2_LIKE / timing-free / ZXZXZX | 0.00099130739442 | 0.0009027778934 | 1.069372 | 0.973870 |
| JONES_MURALI_2026 / SEQ / ZZZXXX | 0.0021817914988 | 0.0018200903519 | 1.071081 | 0.893515 |
| JONES_MURALI_2026 / SEQ / XXXZZZ | 0.0019630253544 | 0.0016379584086 | 1.079178 | 0.900472 |
| JONES_MURALI_2026 / SEQ / ZXZXZX | 0.0024456846756 | 0.0020404169837 | 1.090363 | 0.909682 |
| JONES_MURALI_2026 / PAR_MR / ZZZXXX | 0.0020168313463 | 0.0016910443022 | 1.071073 | 0.898059 |
| JONES_MURALI_2026 / PAR_MR / XXXZZZ | 0.0018683411658 | 0.0015668293773 | 1.056754 | 0.886216 |
| JONES_MURALI_2026 / PAR_MR / ZXZXZX | 0.0022955787993 | 0.0019248823815 | 1.070200 | 0.897381 |

### Exact Uniform(0.001) alternating-order mechanism

All **431280 distinct-location payload pairs per design** are evaluated. Identity matches stage, actual check kind/support, original operation index and payload; execution position is remapped. A3 reproduces the approved raw pair sums exactly within floating-point enumeration precision:

- ZZZXXX P2 = 0.000767911111111; survival-corrected = 0.000698983372468.
- ZXZXZX P2 = 0.000876146666667; survival-corrected = 0.000797503699298.
- Delta P2 = 0.000108235555556; corrected delta = 9.85203268302e-05. Survival factor is (1-0.001)^94 = 0.910239951414, including the approved potential-location convention.

| Matched physical payload-pair partition | Pair count | Raw P2 weight | Corrected weight |
| --- | ---: | ---: | ---: |
| new_only | 18846 | 0.000157182222222 | 0.000143073538319 |
| removed | 9949 | 4.89466666667e-05 | 4.45532114885e-05 |
| common | 82891 | 0.000718964444444 | 0.00065443016098 |

`new_only - removed` equals the net delta; common cases are not double-counted. The complete JSON lists **282 positive classes and 132 negative classes**, each with location/operation types, both circuit-position maps and both trigger paths. This includes countervailing improvements, rather than reporting only newly malignant examples. Pure SPAM pair weights have zero net order penalty here; every positive broad-type excess involves a CX fault.

| Location types | ZZZXXX P2 | ZXZXZX P2 | Delta P2 | Share of net excess |
| --- | ---: | ---: | ---: | ---: |
| 2q*2q | 0.00038884444444 | 0.00042348 | 3.4635555555e-05 | 32.0002% |
| 2q*meas | 0.00017753333333 | 0.00021433333333 | 3.68e-05 | 33.9999% |
| 2q*reset | 0.00017753333333 | 0.00021433333333 | 3.68e-05 | 33.9999% |
| meas*meas | 6e-06 | 6e-06 | 0 | 0.0000% |
| meas*reset | 1.2e-05 | 1.2e-05 | 0 | 0.0000% |
| reset*reset | 6e-06 | 6e-06 | 0 | 0.0000% |

Top positive joint classes (physical identities held fixed, with both actual trigger paths):

| Rank | Types / ops | Physical checks | Positions ZZZXXX → ZXZXZX | Baseline trigger path | Alternating trigger path | Delta P2 | New-only / removed payload cases |
| --- | --- | --- | --- | --- | --- | ---: | ---: |
| 1 | 2q*reset / RX*CX | F:Z2 + B:Z2 | F2,B2 → F4,B4 | F2:Z2:flag1:synd0->bare | F4:Z2:flag1:synd0->bare | 1.4e-06 | 21 / 0 |
| 2 | 2q*meas / MX*CX | F:Z2 + B:Z2 | F2,B2 → F4,B4 | F2:Z2:flag1:synd0->bare | F4:Z2:flag1:synd0->bare | 1.4e-06 | 21 / 0 |
| 3 | 2q*2q / CX*CX | F:X2 + B:Z2 | F5,B2 → F5,B4 | F5:X2:flag0:synd1->bare | F5:X2:flag0:synd1->bare | 1.333333333e-06 | 372 / 72 |
| 4 | 2q*reset / RX*CX | F:Z2 + B:X0 | F2,B3 → F4,B1 | F2:Z2:flag1:synd0->bare | F4:Z2:flag1:synd0->bare | 1.133333333e-06 | 17 / 0 |
| 5 | 2q*meas / MX*CX | F:Z2 + B:X0 | F2,B3 → F4,B1 | F2:Z2:flag1:synd0->bare | F4:Z2:flag1:synd0->bare | 1.133333333e-06 | 17 / 0 |
| 6 | 2q*2q / CX*CX | F:Z2 + B:Z2 | F2,B2 → F4,B4 | F2:Z2:flag0:synd1->bare | F4:Z2:flag0:synd1->bare | 1.12e-06 | 336 / 84 |
| 7 | 2q*reset / RX*CX | F:Z0 + B:X0 | F0,B3 → F0,B1 | F0:Z0:flag1:synd0->bare | F0:Z0:flag1:synd0->bare | 1.066666667e-06 | 16 / 0 |
| 8 | 2q*meas / MX*CX | F:Z0 + B:X0 | F0,B3 → F0,B1 | F0:Z0:flag1:synd0->bare | F0:Z0:flag1:synd0->bare | 1.066666667e-06 | 16 / 0 |
| 9 | 2q*reset / RX*CX | F:Z1 + B:X0 | F1,B3 → F2,B1 | F1:Z1:flag1:synd0->bare | F2:Z1:flag1:synd0->bare | 1.066666667e-06 | 16 / 0 |
| 10 | 2q*meas / MX*CX | F:Z1 + B:X0 | F1,B3 → F2,B1 | F1:Z1:flag1:synd0->bare | F2:Z1:flag1:synd0->bare | 1.066666667e-06 | 16 / 0 |

Top positive position-aligned aggregates (the JSON additionally maps each physical check identity, so equal positions are not mistaken for equal checks):

| Types | Executed positions | ZZZXXX P2 | ZXZXZX P2 | Delta P2 |
| --- | --- | ---: | ---: | ---: |
| 2q*2q | F4 + B3 | 7.653333333e-06 | 1.297777778e-05 | 5.324444444e-06 |
| 2q*reset | F4 + B3 | 4.133333333e-06 | 9.333333333e-06 | 5.2e-06 |
| 2q*meas | F4 + B3 | 4.133333333e-06 | 9.333333333e-06 | 5.2e-06 |
| 2q*2q | F4 + B1 | 7.933333333e-06 | 1.280444444e-05 | 4.871111111e-06 |
| 2q*2q | F1 + B4 | 7.831111111e-06 | 1.268444444e-05 | 4.853333333e-06 |
| 2q*reset | F5 + B2 | 3.6e-06 | 8.333333333e-06 | 4.733333333e-06 |
| 2q*meas | F5 + B2 | 3.6e-06 | 8.333333333e-06 | 4.733333333e-06 |
| 2q*reset | F1 + B4 | 4.333333333e-06 | 8.8e-06 | 4.466666667e-06 |
| 2q*meas | F1 + B4 | 4.333333333e-06 | 8.8e-06 | 4.466666667e-06 |
| 2q*2q | F5 + B2 | 6.915555556e-06 | 1.126666667e-05 | 4.351111111e-06 |

#### Propagation of the top contributors

All qubit and op indices below are **zero-based**; support labels Z2/X0 refer to the original Hamming check rows, not execution positions. Faults act at the existing post-gate/reset locations or flip the reported measurement. Full scalar per-operation frames/readouts, actual executed prefixes, syndromes, corrections and endpoint logical decisions are stored in `top_examples` in `results/stageA3_pair_mechanism.json`.

1. **False flag plus a late bare Z2 hook (top classes 1 and 2).** A Z fault after `RX(8)` at flagged Z2 op 1, or a `FLIP` of its flag `MX(8)` at op 9, produces a spurious flag with syndrome 0. It triggers flagged position F2 in ZZZXXX or F4 in ZXZXZX. Pair it with `IY` after bare Z2 op 2, `CX(2,7)` (B2 versus B4). The syndrome ancilla's Z component propagates through the remaining CNOTs to data **Z4 Z6**; its X component also flips that Z-check readout. Under ZZZXXX all three bare X-checks follow the hook, giving s_z=2; the unchanged flagged table applies Z4 Z6. A residual X0 from s_x=1 is weight one and is removed by the final ideal decoder. Under ZXZXZX, X0 and X1 were read before the hook and X2 sees even parity, giving stale/incomplete s_z=0 and no Z correction. The final residual contains Z4 Z6; ideal decoding adds Z1 and leaves logical Z1 Z4 Z6. The reset and measurement classes each add 21 newly malignant payload combinations with no removals, delta P2=1.4e-6 per class.

2. **False flag plus an X4 error between bare Z readouts (classes 4 and 5).** Use the same spurious flagged Z2 signal, paired with `IX` after bare X0 op 2, `CX(7,4)` (B3 in grouped order, B1 in alternating order). It leaves data X4. Grouped ZZZXXX has already read all bare Z checks, so s_x=0 and X4 remains a correctable single error at the endpoint. Alternating ZXZXZX has read Z0 but not Z1/Z2; subsequent readouts produce incomplete s_x=1. The standard X decoder applies X0, converting X4 into X0 X4, and the final ideal decoder adds X3, leaving logical X0 X3 X4. This example directly shows the one-round boundary asymmetry: one design defers a weight-one error, while the other creates an inconsistent correction. It is not a single-fault FT violation, since the false trigger and gate fault are two faults. Each corresponding reset/measurement class contributes 17 new payload cases with no removals, delta P2=1.133333e-6.

3. **Two-CX stale-syndrome case (class 3).** A `YI` fault after flagged X2 op 2, `CX(7,0)`, is paired with `YI` after bare Z2 op 2, `CX(2,7)`. The first makes the actual F5 X2 syndrome trigger (flag 0); the second adds a late Z2 component during the bare round. Both designs have s_x=1, but grouped order reads all bare X checks afterward and obtains s_z=3, applying Z2. Alternating order has already read X0/X1 and obtains s_z=1, applying Z0 instead. Its Z0 Z2 residual becomes logical Z0 Z1 Z2 after ideal decoding. This joint class contains 372 newly malignant and 72 removed payload combinations, net delta P2=1.333333e-6; the removals are retained in the class accounting.

These traces explain the main mechanism as **syndrome bits collected on different sides of a later data error**, together with the unchanged flagged-table/standard-decoder response. It is not an invented extra idle penalty: this attribution uses timing-free uniform noise and identical physical fault pairs. Carrying errors into later rounds removes the free boundary treatment for some cases, but does not imply that all alternating-order losses or adaptive ordering differences must vanish.

For the one-round mechanism attribution, approved high-shot one-round MC supplies the reference; its raw/corrected P2 comparisons are stored in the mechanism JSON. Epsilon_tot remains primary for design comparisons. The full noisy memory MC includes cross-round faults and higher-order events.

Assumptions made:
- The user's A3 request overrides the earlier one-round-only scope. No added waiting, hardware timing, layout, recooling or transport is introduced; H2 retains timing-free channels.
- Actual protocol corrections occur each round. The frame is carried throughout and ideal decoding is only at the last endpoint. Logical changes can cancel; only an explicit preparation abort would be absorbing for a future adapter.
- Per the user's clarification, primary fits are epsilon_X/Z marginals, and every primary design contrast uses the derived epsilon_tot with joint propagated uncertainty. Union data remain descriptive, including R=1; no union decay is fitted.
- R2 error bars mean one-sigma MC stderr, with separate combined trend uncertainty shown. All R≥4 fits are supplied, not selected only when they improve conclusions. Poor goodness of fit limits stationary interpretation even if statistical precision is high.
- Three-sigma statements concern the specified independent-design sampling contrast. Fit quality and heuristic dispersion diagnostics do not authorize an infinite-R conclusion or a bias correction.
- One-round P2 is always accompanied by exact survival correction and comparison to one-round MC. No multi-round P2×R approximation or repurposed one-round memory p_eff is claimed.
- Stage B3's primary metric is epsilon_tot from these marginal memory fits; one-round union p_L is secondary. The reusable harness is ready, but no Steane-style ancilla/circuit work has started.

Open questions for the human:
- Review the persistent epsilon_Z difference and the smaller epsilon_tot difference for grouped orders, including the explicitly reported R≥4 results and fit quality. The tested finite windows do not establish grouped-order stationary equivalence.
- Some fits remain poor after excluding R=2; the report preserves this limitation rather than silently treating a precise finite-window estimate as an asymptotic rate.
- Stop here. Stage B1 and all later Stage B work await approval.


## Stage B1 — Published verified ancilla preparation

Completed B1 only. No full Steane EC round, data/ancilla coupling, B2 decoder, or B3 design comparison was assembled or sampled. The later B3 decisions are recorded below.

### Circuit, derivation and single-fault checks

The fixed gadget is Goto, Scientific Reports **6**, 19578 (2016), [Fig. 1(c)](https://doi.org/10.1038/srep19578); the prompt's arXiv:1610.03309 identifies an unrelated paper. The [publisher correction](https://doi.org/10.1038/s41598-022-06633-6) corrects the volume. [Peham et al., Example 3/Fig. 4](https://doi.org/10.1103/PRXQuantum.6.020330) supplies a further published circuit/dual reference. The annotated implementation is `figures/verified_prep_circuits.png`.

Block-local qubits 0..6 map to physical 7..13; verifier is 14, with data 0..6 reserved and untouched. Start the zero gadget with plus inputs at local 1,2,3 and zero inputs elsewhere. Encoder CNOTs, in execution order, are (1,0), (3,5), (2,6), (1,4), (2,0), (3,6), (1,5), (6,4). The generated X checks have supports 0145, 0246, 3456: their row span equals the frozen code. Verify logical Z on local 0,5,6 using three CNOTs into the verifier and Z readout; accept reported 0. The plus gadget is its complete Hadamard dual: swap reset/readout bases and reverse every CNOT, without adding gates. Each attempt has eight resets, eleven CNOTs and one readout.

Among 141 encoder single-fault payloads, 30 have dangerous weight above one before verification. Modulo code stabilizers these form two cosets, represented by X0X2 (22 payloads) and X0X4 (8); both anticommute with the verifier. The plus case is dual. Dangerous X for zero / Z for plus is reduced modulo code stabilizers only, never modulo a harmful logical operator. All accepted dangerous components have reduced weight at most one.

Universal checks: **1302 faults across six basis/timing geometries; zero violations and zero frame/tableau mismatches**. Positive-channel checks: **3288 faults across all 18 rows; zero violations and mismatches**. The maximum accepted orthogonal component has raw and code-reduced weight **3**, but prepared-state-reduced weight **1** (logical Z stabilizes zero; logical X stabilizes plus). Orthogonal Z on the zero control or X on the plus target does not propagate to the future data block through the intended CNOT direction. This is a component-level preparation result; B2 must still verify the complete EC gadget. Exact failing-fault arrays are empty; no circuit or decoder tuning followed these checks.

### Per-attempt Monte Carlo and pair attribution

Each row uses **5,000,000 independent attempts**. Acceptance is the primary B1 metric. The hazard below is the unconditional probability of verifier acceptance with dangerous component reduced weight >1; it is **not full-EC p_L or memory epsilon**. Raw P2, exact survival-corrected P2 and the mean-rate correction are attribution diagnostics for that same hazard. The survival factor is the product of (1-p_loc) over all other locations. The JSON additionally reports conditional hazards, Wilson 95% intervals, exact finite-rate convolution and all P2/MC ratios and differences in MC standard errors.

| Preset / policy | State | Acceptance MC ± 1σ | Hazard MC ± 1σ | Raw P2 | Survival P2 | Mean survival P2 | Raw / MC | Survival / MC |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| UniformNoise(0.0005) / timing-free | 0 | 0.99459980 ± 3.3e-05 | 6.2e-06 ± 1.1e-06 | 5.181111e-06 | 5.134679e-06 | 5.134679e-06 | 0.836 | 0.828 |
| UniformNoise(0.0005) / timing-free | + | 0.99461600 ± 3.3e-05 | 5e-06 ± 1e-06 | 5.181111e-06 | 5.134679e-06 | 5.134679e-06 | 1.036 | 1.027 |
| UniformNoise(0.001) / timing-free | 0 | 0.98921420 ± 4.6e-05 | 1.88e-05 ± 1.9e-06 | 2.072444e-05 | 2.035456e-05 | 2.035456e-05 | 1.102 | 1.083 |
| UniformNoise(0.001) / timing-free | + | 0.98922700 ± 4.6e-05 | 1.98e-05 ± 2e-06 | 2.072444e-05 | 2.035456e-05 | 2.035456e-05 | 1.047 | 1.028 |
| UniformNoise(0.0015) / timing-free | 0 | 0.98385520 ± 5.6e-05 | 5.04e-05 ± 3.2e-06 | 4.663e-05 | 4.538691e-05 | 4.538691e-05 | 0.925 | 0.901 |
| UniformNoise(0.0015) / timing-free | + | 0.98398860 ± 5.6e-05 | 5.14e-05 ± 3.2e-06 | 4.663e-05 | 4.538691e-05 | 4.538691e-05 | 0.907 | 0.883 |
| TWO_QUBIT_ONLY(0.0005) / timing-free | 0 | 0.99707240 ± 2.4e-05 | 1.8e-06 ± 6e-07 | 1.831111e-06 | 1.822888e-06 | 1.822068e-06 | 1.017 | 1.013 |
| TWO_QUBIT_ONLY(0.0005) / timing-free | + | 0.99705940 ± 2.4e-05 | 3e-06 ± 7.7e-07 | 1.831111e-06 | 1.822888e-06 | 1.822068e-06 | 0.610 | 0.608 |
| TWO_QUBIT_ONLY(0.001) / timing-free | 0 | 0.99416980 ± 3.4e-05 | 7.4e-06 ± 1.2e-06 | 7.324444e-06 | 7.258788e-06 | 7.25227e-06 | 0.990 | 0.981 |
| TWO_QUBIT_ONLY(0.001) / timing-free | + | 0.99421600 ± 3.4e-05 | 6e-06 ± 1.1e-06 | 7.324444e-06 | 7.258788e-06 | 7.25227e-06 | 1.221 | 1.210 |
| TWO_QUBIT_ONLY(0.002) / timing-free | 0 | 0.98848060 ± 4.8e-05 | 3.42e-05 ± 2.6e-06 | 2.929778e-05 | 2.877462e-05 | 2.872307e-05 | 0.857 | 0.841 |
| TWO_QUBIT_ONLY(0.002) / timing-free | + | 0.98841720 ± 4.8e-05 | 3.02e-05 ± 2.5e-06 | 2.929778e-05 | 2.877462e-05 | 2.872307e-05 | 0.970 | 0.953 |
| H2_LIKE / timing-free | 0 | 0.99142240 ± 4.1e-05 | 1.04e-05 ± 1.4e-06 | 1.486981e-05 | 1.467745e-05 | 1.466391e-05 | 1.430 | 1.411 |
| H2_LIKE / timing-free | + | 0.99129660 ± 4.2e-05 | 1.1e-05 ± 1.5e-06 | 1.515138e-05 | 1.4954e-05 | 1.494054e-05 | 1.377 | 1.359 |
| JONES_MURALI_2026 / SEQ | 0 | 0.97362660 ± 7.2e-05 | 0.0001122 ± 4.7e-06 | 0.0001068578 | 0.0001020384 | 0.0001014549 | 0.952 | 0.909 |
| JONES_MURALI_2026 / SEQ | + | 0.97311880 ± 7.2e-05 | 0.0001252 ± 5e-06 | 0.0001160877 | 0.0001108149 | 0.0001102124 | 0.927 | 0.885 |
| JONES_MURALI_2026 / PAR_MR | 0 | 0.97373420 ± 7.2e-05 | 0.0001066 ± 4.6e-06 | 0.0001068578 | 0.000102071 | 0.000101492 | 1.002 | 0.958 |
| JONES_MURALI_2026 / PAR_MR | + | 0.97331360 ± 7.2e-05 | 0.0001122 ± 4.7e-06 | 0.0001151445 | 0.0001099553 | 0.0001093659 | 1.026 | 0.980 |

No acceptance test uses ideal error weight. All multi-fault verifier-zero states are accepted, including dangerous states. Exact XOR convolution keeps all fault orders and is a secondary benchmark. Extra orthogonal idle faults can leave a hazard unchanged while excluding that history from exactly-two-fault survival P2; corrected P2 therefore need not approximate the full hazard better.

The original H2 rare-hazard MC samples (52 and 55 events) are low compared with exact expectations (73.85 and 75.24). They remain the primary recorded samples. Fixed independent **20M-attempt held-out checks per basis** gave 326 and 300 events versus expectations 295.40 and 300.96 (null-variance z=+1.78 and −0.056). These are reported separately, not substituted for the initial samples. All scalar/frame/tableau checks pass.

### Verification branch triggers

The verification rejection probability triggers a fresh preparation attempt. Disjoint per-attempt rates below classify the nonzero parity contributions by type; mixed means more than one type contributes. Individual contributing-type and fractional allocations are also in JSON. Flag-trigger rate and failure share through a bare round are **not applicable** to B1; no flag or bare round runs here. Bare-only versus flagged single-fault counts are likewise inapplicable: the table above covers preparation faults, with reset/encoder/verification/readout counts retained per row.

| Preset / policy | State | Rejection | Reset only | Measurement only | 2q only | Idle only | Mixed |
|---|---|---:|---:|---:|---:|---:|---:|
| UniformNoise(0.0005) / timing-free | 0 | 0.0054002 | 0.0019528 | 0.000493 | 0.0029544 | 0 | 0 |
| UniformNoise(0.0005) / timing-free | + | 0.005384 | 0.0019668 | 0.0005014 | 0.0029158 | 0 | 0 |
| UniformNoise(0.001) / timing-free | 0 | 0.0107858 | 0.0039672 | 0.0009758 | 0.0058428 | 0 | 0 |
| UniformNoise(0.001) / timing-free | + | 0.010773 | 0.0039742 | 0.000974 | 0.0058248 | 0 | 0 |
| UniformNoise(0.0015) / timing-free | 0 | 0.0161448 | 0.0059678 | 0.0014918 | 0.0086848 | 0 | 4e-07 |
| UniformNoise(0.0015) / timing-free | + | 0.0160114 | 0.0058802 | 0.0014766 | 0.0086544 | 0 | 2e-07 |
| TWO_QUBIT_ONLY(0.0005) / timing-free | 0 | 0.0029276 | 0 | 0 | 0.0029276 | 0 | 0 |
| TWO_QUBIT_ONLY(0.0005) / timing-free | + | 0.0029406 | 0 | 0 | 0.0029406 | 0 | 0 |
| TWO_QUBIT_ONLY(0.001) / timing-free | 0 | 0.0058302 | 0 | 0 | 0.0058302 | 0 | 0 |
| TWO_QUBIT_ONLY(0.001) / timing-free | + | 0.005784 | 0 | 0 | 0.005784 | 0 | 0 |
| TWO_QUBIT_ONLY(0.002) / timing-free | 0 | 0.0115194 | 0 | 0 | 0.0115194 | 0 | 0 |
| TWO_QUBIT_ONLY(0.002) / timing-free | + | 0.0115828 | 0 | 0 | 0.0115828 | 0 | 0 |
| H2_LIKE / timing-free | 0 | 0.0085776 | 0.0001514 | 0.0009868 | 0.0074394 | 0 | 0 |
| H2_LIKE / timing-free | + | 0.0087034 | 0.0002276 | 0.0010248 | 0.007451 | 0 | 0 |
| JONES_MURALI_2026 / SEQ | 0 | 0.0263734 | 0.0196462 | 0.0009956 | 0.0057312 | 0 | 4e-07 |
| JONES_MURALI_2026 / SEQ | + | 0.0268812 | 0.0196916 | 0.001015 | 0.0057128 | 0.0004618 | 0 |
| JONES_MURALI_2026 / PAR_MR | 0 | 0.0262658 | 0.0195798 | 0.0009908 | 0.0056948 | 0 | 4e-07 |
| JONES_MURALI_2026 / PAR_MR | + | 0.0266864 | 0.0196848 | 0.000965 | 0.0057232 | 0.000313 | 4e-07 |

### Timing and bounded retry diagnostics

Only sourced Jones–Murali timing is used. H2 remains timing-free (Stage 1 channels only); transport is empty and no hardware schedule is invented. The eight reset operations overlap in PAR_MR; encoder CNOTs remain serial, and there is just one verifier readout. A block ion becomes live after its reset finishes and waits through verifier readout. Single-attempt durations (zero, plus): SEQ **1.255, 1.270 ms**; PAR_MR **0.895, 0.900 ms**. Per-ion idle seconds and summed Z probabilities are recorded in each timed JSON row.

These are standalone eight-ion preparations in a reserved fifteen-ion layout. **External data waiting is outside B1 and must be added in B2 on every attempt, including rejected attempts.** No full-EC time or resource claim is made here. Rejected attempts are discarded and all eight ions reset. Accepted states retain full Pauli frames for a future B2 adapter.

The user resolved the retry policy as **max_attempts = 4**, an initial attempt plus three retries. The K=4 result is nominal; K=3 is retained as a standalone diagnostic. Each cap/row has 5M independently sampled cycles. For rejection q, exact exhaustion is q^K and mean attempts is Σ(q^i), i=0..K−1. Exhaustion counts as a standalone failure; the JSON also includes accepted dangerous states, total standalone failures and resource counts including rejected attempts.

| Preset / policy | State | K | Exhaustion MC | Exact exhaustion | Mean attempts MC | Mean elapsed ms (timed only) |
|---|---|---:|---:|---:|---:|---:|
| UniformNoise(0.0005) / timing-free | 0 | 3 | 3/5000000 | 1.579833e-07 | 1.0054604 | — |
| UniformNoise(0.0005) / timing-free | 0 | 4 | 0/5000000 | 8.540469e-10 | 1.0054142 | — |
| UniformNoise(0.0005) / timing-free | + | 3 | 0/5000000 | 1.579833e-07 | 1.0054126 | — |
| UniformNoise(0.0005) / timing-free | + | 4 | 0/5000000 | 8.540469e-10 | 1.0054592 | — |
| UniformNoise(0.001) / timing-free | 0 | 3 | 6/5000000 | 1.244863e-06 | 1.0108800 | — |
| UniformNoise(0.001) / timing-free | 0 | 4 | 0/5000000 | 1.339148e-08 | 1.0108280 | — |
| UniformNoise(0.001) / timing-free | + | 3 | 5/5000000 | 1.244863e-06 | 1.0108938 | — |
| UniformNoise(0.001) / timing-free | + | 4 | 0/5000000 | 1.339148e-08 | 1.0109156 | — |
| UniformNoise(0.0015) / timing-free | 0 | 3 | 24/5000000 | 4.138313e-06 | 1.0162322 | — |
| UniformNoise(0.0015) / timing-free | 0 | 4 | 0/5000000 | 6.644022e-08 | 1.0162502 | — |
| UniformNoise(0.0015) / timing-free | + | 3 | 16/5000000 | 4.138313e-06 | 1.0162788 | — |
| UniformNoise(0.0015) / timing-free | + | 4 | 0/5000000 | 6.644022e-08 | 1.0162950 | — |
| TWO_QUBIT_ONLY(0.0005) / timing-free | 0 | 3 | 0/5000000 | 2.503865e-08 | 1.0029018 | — |
| TWO_QUBIT_ONLY(0.0005) / timing-free | 0 | 4 | 0/5000000 | 7.325115e-11 | 1.0029182 | — |
| TWO_QUBIT_ONLY(0.0005) / timing-free | + | 3 | 1/5000000 | 2.503865e-08 | 1.0029108 | — |
| TWO_QUBIT_ONLY(0.0005) / timing-free | + | 4 | 0/5000000 | 7.325115e-11 | 1.0029506 | — |
| TWO_QUBIT_ONLY(0.001) / timing-free | 0 | 3 | 0/5000000 | 1.987144e-07 | 1.0057488 | — |
| TWO_QUBIT_ONLY(0.001) / timing-free | 0 | 4 | 0/5000000 | 1.159593e-09 | 1.0059334 | — |
| TWO_QUBIT_ONLY(0.001) / timing-free | + | 3 | 0/5000000 | 1.987144e-07 | 1.0059002 | — |
| TWO_QUBIT_ONLY(0.001) / timing-free | + | 4 | 0/5000000 | 1.159593e-09 | 1.0058438 | — |
| TWO_QUBIT_ONLY(0.002) / timing-free | 0 | 3 | 11/5000000 | 1.564522e-06 | 1.0117466 | — |
| TWO_QUBIT_ONLY(0.002) / timing-free | 0 | 4 | 0/5000000 | 1.81625e-08 | 1.0117306 | — |
| TWO_QUBIT_ONLY(0.002) / timing-free | + | 3 | 5/5000000 | 1.564522e-06 | 1.0116394 | — |
| TWO_QUBIT_ONLY(0.002) / timing-free | + | 4 | 0/5000000 | 1.81625e-08 | 1.0116778 | — |
| H2_LIKE / timing-free | 0 | 3 | 7/5000000 | 6.361984e-07 | 1.0086672 | — |
| H2_LIKE / timing-free | 0 | 4 | 0/5000000 | 5.471714e-09 | 1.0087072 | — |
| H2_LIKE / timing-free | + | 3 | 3/5000000 | 6.575131e-07 | 1.0087384 | — |
| H2_LIKE / timing-free | + | 4 | 0/5000000 | 5.717496e-09 | 1.0087974 | — |
| JONES_MURALI_2026 / SEQ | 0 | 3 | 79/5000000 | 1.810219e-05 | 1.0269024 | 1.288763 |
| JONES_MURALI_2026 / SEQ | 0 | 4 | 2/5000000 | 4.753077e-07 | 1.0269538 | 1.288827 |
| JONES_MURALI_2026 / SEQ | + | 3 | 90/5000000 | 1.90632e-05 | 1.0275094 | 1.304937 |
| JONES_MURALI_2026 / SEQ | + | 4 | 2/5000000 | 5.092462e-07 | 1.0275590 | 1.305000 |
| JONES_MURALI_2026 / PAR_MR | 0 | 3 | 89/5000000 | 1.810219e-05 | 1.0270486 | 0.919208 |
| JONES_MURALI_2026 / PAR_MR | 0 | 4 | 2/5000000 | 4.753077e-07 | 1.0269586 | 0.919128 |
| JONES_MURALI_2026 / PAR_MR | + | 3 | 86/5000000 | 1.881783e-05 | 1.0272222 | 0.924500 |
| JONES_MURALI_2026 / PAR_MR | + | 4 | 3/5000000 | 5.005254e-07 | 1.0273836 | 0.924645 |

Zero observed exhaustions do not imply zero probability; Wilson intervals and nonzero exact expectations remain in JSON. Neither retry convention changes the per-attempt acceptance results.

### Files and validation

Added `src/verified_prep.py` (published circuits, independent tableau/frame checks, finite-rate convolution, pair counting and retry/frame APIs), `src/prep_experiment.py` (resumable experiments/report), `tests/test_verified_prep.py`, and the A3 approved-artifact hash fixture. Results are `results/stageB1_preparation.json`, `results/stageB1_validation.json`, and resumable `results/stageB1_cache/`; diagram is `figures/verified_prep_circuits.png`. README and noise-model documentation describe reproduction, timing boundaries and scope.

Validation: **143 tests passed in 414.86 s**, including all prior regression tests and 24 new B1 checks. All **147 approved artifacts are byte-identical**. This preserves the approved single-fault counts, pair counts, seeded MC, A-stage fits/results and synthetic transport outputs. Recorded primary sampling is 90M attempts, retry sampling is 180,000,000 cycles / 182,323,349 actual attempts, plus 40M held-out H2 attempts. Count conservation, cause allocations, pair attribution sums, survival bounds and retry identities pass; both circuit panels were visually inspected. No dependencies were added.


FT checks and cross-check counts (positive channel payloads):

| Preset / policy | State | Fault runs | Accepted | Rejected | FT failures | Frame/tableau mismatches |
|---|---|---:|---:|---:|---:|---:|
| UniformNoise(0.0005) / timing-free | 0 | 174 | 81 | 93 | 0 | 0 |
| UniformNoise(0.0005) / timing-free | + | 174 | 81 | 93 | 0 | 0 |
| UniformNoise(0.001) / timing-free | 0 | 174 | 81 | 93 | 0 | 0 |
| UniformNoise(0.001) / timing-free | + | 174 | 81 | 93 | 0 | 0 |
| UniformNoise(0.0015) / timing-free | 0 | 174 | 81 | 93 | 0 | 0 |
| UniformNoise(0.0015) / timing-free | + | 174 | 81 | 93 | 0 | 0 |
| TWO_QUBIT_ONLY(0.0005) / timing-free | 0 | 165 | 77 | 88 | 0 | 0 |
| TWO_QUBIT_ONLY(0.0005) / timing-free | + | 165 | 77 | 88 | 0 | 0 |
| TWO_QUBIT_ONLY(0.001) / timing-free | 0 | 165 | 77 | 88 | 0 | 0 |
| TWO_QUBIT_ONLY(0.001) / timing-free | + | 165 | 77 | 88 | 0 | 0 |
| TWO_QUBIT_ONLY(0.002) / timing-free | 0 | 165 | 77 | 88 | 0 | 0 |
| TWO_QUBIT_ONLY(0.002) / timing-free | + | 165 | 77 | 88 | 0 | 0 |
| H2_LIKE / timing-free | 0 | 180 | 87 | 93 | 0 | 0 |
| H2_LIKE / timing-free | + | 184 | 87 | 97 | 0 | 0 |
| JONES_MURALI_2026 / SEQ | 0 | 222 | 129 | 93 | 0 | 0 |
| JONES_MURALI_2026 / SEQ | + | 226 | 111 | 115 | 0 | 0 |
| JONES_MURALI_2026 / PAR_MR | 0 | 220 | 127 | 93 | 0 | 0 |
| JONES_MURALI_2026 / PAR_MR | + | 222 | 110 | 112 | 0 | 0 |

Frozen-result checks (all pass in the full suite, with approved artifact hashes unchanged):

| Regression | Result |
|---|---|
| Flag single faults | 0/1218 |
| No-flag single faults | 48/1218 |
| Cross-check faults | 0/1176 mismatches |
| Malignant pairs / A dictionaries | 92840 / exact match |
| Seeded Uniform MC | 542/3M, seed 2026 |
| Seeded Jones–Murali MC | 6121/3M, seed 3326 |
| Stage 4 empty transport and synthetic monotonicity | Pass |
| A1/A2/A3 fixtures and frozen outputs | Pass; 147 artifact hashes unchanged |

Assumptions made:

- The published zero circuit and full Hadamard dual are fixed; paper-qubit mapping is explicit and the frozen stabilizer row span is preserved.
- Reset basis durations and idle channels reuse the approved model; data storage during preparation is deferred to B2, where it belongs in the coupled gadget. H2 receives no invented timing.
- Preparation-hazard P2 is not full-EC p_L. B1 does not produce p_eff or memory epsilon. Future B3 primary epsilon_tot comes from the A3 memory harness, with epsilon_X and epsilon_Z always shown.
- Later B3 uses XXXZZZ as the flag best-order arm under every preset, including nonsignificant cases, and both Steane half-orders (X-syndrome first and Z-syndrome first). These are recorded decisions only; no B3 action was taken.

Open questions for the human:

- Review the B1 circuit, acceptance and component-level fault-tolerance results.
- Retry policy resolved by the user: `max_attempts = 4`. The clarification table below reports r and r^4.
- **Stop here. B2 and all later work require the next checkpoint approval.**


### B1 clarification — max_attempts = 4

The user specifies an initial attempt plus three retries, **four attempts total**. The runner uses `max_attempts = 4`; the CLI is `--max-attempts 4`. Previously sampled three-attempt diagnostics remain available. All original MC samples, pair results, single-fault checks and caches are unchanged. No B2/B3 execution is authorized by this clarification.

For independent fresh attempts, per-attempt rejection r implies all four attempts fail with probability r^4 for **one standalone preparation**. MC-derived r and r^4 below are primary; uncertainty on r^4 is propagated as 4r³σ_r, and transformed Wilson 95% bounds are in JSON. The exact finite-rate r^4 is a secondary benchmark. This is not the abort rate of a complete two-half EC round or a multi-round memory experiment.

| Preset / policy | State | Rejection r MC ± 1σ | Implied r^4 ± 1σ | Exact r^4 | Observed four-attempt exhaustion |
|---|---|---:|---:|---:|---:|
| UniformNoise(0.0005) / timing-free | 0 | 0.00540020 ± 3.3e-05 | 8.5043158e-10 ± 2.1e-11 | 8.5404694e-10 | 0/5000000 |
| UniformNoise(0.0005) / timing-free | + | 0.00538400 ± 3.3e-05 | 8.4027261e-10 ± 2e-11 | 8.5404694e-10 | 0/5000000 |
| UniformNoise(0.001) / timing-free | 0 | 0.01078580 ± 4.6e-05 | 1.3533479e-08 ± 2.3e-10 | 1.3391481e-08 | 0/5000000 |
| UniformNoise(0.001) / timing-free | + | 0.01077300 ± 4.6e-05 | 1.346935e-08 ± 2.3e-10 | 1.3391481e-08 | 0/5000000 |
| UniformNoise(0.0015) / timing-free | 0 | 0.01614480 ± 5.6e-05 | 6.7940803e-08 ± 9.5e-10 | 6.644022e-08 | 0/5000000 |
| UniformNoise(0.0015) / timing-free | + | 0.01601140 ± 5.6e-05 | 6.5722977e-08 ± 9.2e-10 | 6.644022e-08 | 0/5000000 |
| TWO_QUBIT_ONLY(0.0005) / timing-free | 0 | 0.00292760 ± 2.4e-05 | 7.3459328e-11 ± 2.4e-12 | 7.3251147e-11 | 0/5000000 |
| TWO_QUBIT_ONLY(0.0005) / timing-free | + | 0.00294060 ± 2.4e-05 | 7.4772829e-11 ± 2.5e-12 | 7.3251147e-11 | 0/5000000 |
| TWO_QUBIT_ONLY(0.001) / timing-free | 0 | 0.00583020 ± 3.4e-05 | 1.1554039e-09 ± 2.7e-11 | 1.1595932e-09 | 0/5000000 |
| TWO_QUBIT_ONLY(0.001) / timing-free | + | 0.00578400 ± 3.4e-05 | 1.119214e-09 ± 2.6e-11 | 1.1595932e-09 | 0/5000000 |
| TWO_QUBIT_ONLY(0.002) / timing-free | 0 | 0.01151940 ± 4.8e-05 | 1.7608381e-08 ± 2.9e-10 | 1.8162497e-08 | 0/5000000 |
| TWO_QUBIT_ONLY(0.002) / timing-free | + | 0.01158280 ± 4.8e-05 | 1.7999243e-08 ± 3e-10 | 1.8162497e-08 | 0/5000000 |
| H2_LIKE / timing-free | 0 | 0.00857760 ± 4.1e-05 | 5.4133133e-09 ± 1e-10 | 5.4717144e-09 | 0/5000000 |
| H2_LIKE / timing-free | + | 0.00870340 ± 4.2e-05 | 5.737937e-09 ± 1.1e-10 | 5.7174958e-09 | 0/5000000 |
| JONES_MURALI_2026 / SEQ | 0 | 0.02637340 ± 7.2e-05 | 4.8379847e-07 ± 5.3e-09 | 4.7530775e-07 | 2/5000000 |
| JONES_MURALI_2026 / SEQ | + | 0.02688120 ± 7.2e-05 | 5.2214919e-07 ± 5.6e-09 | 5.0924623e-07 | 2/5000000 |
| JONES_MURALI_2026 / PAR_MR | 0 | 0.02626580 ± 7.2e-05 | 4.7595132e-07 ± 5.2e-09 | 4.7530775e-07 | 2/5000000 |
| JONES_MURALI_2026 / PAR_MR | + | 0.02668640 ± 7.2e-05 | 5.0717748e-07 ± 5.5e-09 | 5.0052544e-07 | 3/5000000 |

Recorded for B2/B3, with no action now:

- Aborts count as logical failures. Report the observed abort rate and its share of epsilon_tot; epsilon_X and epsilon_Z remain mandatory.
- Run Jones–Murali SEQ cap sensitivity at max_attempts ∈ {2, 4, ∞}; the unbounded case retries until acceptance.
- Every rejected attempt contributes real time, fault locations and data waiting.
- Preserve the approved A3 memory metric, XXXZZZ arm under every preset, and both Steane half-orders.

Validation: all 18 r^4 entries and exact cap-four identities checked; all 75 MC/check/analysis caches are byte-identical; all 147 pre-B1 approved artifacts remain byte-identical. The focused B1 suite passes: **24 tests in 3.37 s**; `git diff --check` passes. This changes reporting and the experiment-runner default only. **Stop before B2.**
