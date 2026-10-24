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
  ion_noise.py              Stages 1–2: channels, checks, Monte Carlo and trigger metrics
  ion_timing.py             Sequential timing and lumped idle Z channels
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

1. First p_eff from Monte Carlo, with P2 for a secondary check.
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
# Default: Jones–Murali Stage 1 reference + Stage 2 + T2 sweep,
# and H2-like Stage 1 only. Monte Carlo is the primary p_L estimate.
.venv/bin/python src/ion_noise.py
.venv/bin/python src/ion_noise.py --preset JONES_MURALI_2026
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

Stage 3 will define p_eff primarily from Monte Carlo; the P2-based p_eff will
be a secondary check. It has not been implemented at Stage 2.

Nominal Jones–Murali Stage 2 result (T2 = 2.2 s):

- Monte Carlo p_L = 0.002007 ± 2.58e-05 (3M shots).
- Trigger rate = 10.50293%; 98.7211% of failures went through the bare round.
- Attribution: raw P2 = 0.0021817915, survival-corrected P2 = 0.0018200904.
- Flagged single faults: 0/1332; bare-only (skipped on noiseless input): 0/888; incoming errors: 0/42.
- Full metrics, trigger category splits, and T2 sweep: [Stage 2 report](docs/STAGE_LOG.md#stage-2--sequential-idle-z-noise-and-primary-monte-carlo-reporting).
