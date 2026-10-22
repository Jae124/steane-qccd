# steane-qccd

Steane [[7,1,3]] fault-tolerant gadgets, their logical error rates, and (next) their
schedules on a QCCD trapped-ion architecture. Master's thesis work, Kubiatowicz group,
UC Berkeley.

## Layout

```
src/
  steane_ft.py              Layer 1: flagged EC round (Chao–Reichardt) + single-fault test
  steane_noise.py           Layer 2: per-location noise model, fault-pair counting, Monte Carlo
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

Idle error, shuttling / split / merge error, noisy state preparation. Noise is the same p
for every operation. Pseudothreshold not yet checked against Chao–Reichardt.

## Next

1. Idle locations with a timing model; published trapped-ion rates; first p_eff.
2. Other gadgets: Goto-verified |0_L⟩, transversal H / S / CNOT, logical measurement.
3. QCCD schedule per gadget → idle and shuttle error per location.
