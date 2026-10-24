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
