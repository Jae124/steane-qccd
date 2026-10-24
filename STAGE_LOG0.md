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
