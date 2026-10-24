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
  ion_noise.py              Stage 1: ion channels, checks, pair/Monte Carlo report
  plot_ec_threshold.py      p_L vs p figure
  draw_flagged_circuit.py   Clean drawing of the flagged Z4Z5Z6Z7 circuit
results/                    Generated data (.json, .stim)
figures/                    Generated figures (.png, .svg)
```

## Run

All scripts run from `src/` and write into `results/` and `figures/`.

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

Idle error, shuttling / split / merge error, and noisy logical encoding. Uniform
noise remains the default; operation-specific ion noise is opt-in. Pseudothreshold not yet checked against Chao–Reichardt.

## Next

1. Idle locations with a timing model, then first p_eff.
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
.venv/bin/python src/ion_noise.py
# Or select one preset (also works as `python ion_noise.py` from src/):
.venv/bin/python src/ion_noise.py --preset H2_LIKE --output results/h2_stage1.json
.venv/bin/python src/ion_noise.py --preset JONES_MURALI_2026 --output results/jones_stage1.json
```

The CLI verifies channel fault tolerance and frame/tableau agreement before
counting pairs and sampling. It writes `results/stage1_ion_noise.json` by default.
At Stage 1 it models operation errors only; idle/transport and p_eff are later
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
