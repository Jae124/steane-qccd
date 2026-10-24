"""Centralized physical inputs. Durations are in seconds; probabilities in [0, 1]."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class IonParams:
    name: str
    p_2q: float
    p_1q: float
    p_reset: float
    p_meas: float
    T2: float | None
    ms_duration: float | None
    rotation_duration: float | None
    measurement_duration: float | None
    reset_duration: float | None
    shuttle_duration: float | None
    split_duration: float | None
    merge_duration: float | None
    junction_duration: float | None
    source: str
    caveat: str

    def __post_init__(self):
        for field in ('p_2q', 'p_1q', 'p_reset', 'p_meas'):
            value = getattr(self, field)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f'{field} must be a probability')


# Jones & Murali, arXiv:2510.23519v1, Table 1 and Section 5.1.
# Keep the requested preset name; the arXiv v1 publication date is October 2025.
JONES_MURALI_2026 = IonParams(
    name='JONES_MURALI_2026',
    p_2q=1e-3,  # PLACEHOLDER requested for this study; not their heating model.
    p_1q=1e-5,  # PLACEHOLDER requested for this study; not their heating model.
    p_reset=5e-3,  # Section 5.1, reset bit-flip probability.
    p_meas=1e-3,  # Section 5.1, measurement error probability.
    T2=2.2,  # Section 5.1, coherence time; used from Stage 2 onward.
    ms_duration=40e-6,  # Table 1, two-qubit MS gate.
    rotation_duration=5e-6,  # Table 1, one-qubit rotation.
    measurement_duration=400e-6,  # Table 1, measurement.
    reset_duration=50e-6,  # Table 1, reset.
    shuttle_duration=5e-6,  # Table 1, shuttle.
    split_duration=80e-6,  # Table 1, split.
    merge_duration=80e-6,  # Table 1, merge.
    junction_duration=100e-6,  # Table 1, junction entry/exit.
    source='https://arxiv.org/html/2510.23519v1#S5',
    caveat='Gate probabilities are placeholders pending a heating model.',
)

# LeBlond et al., arXiv:2508.14227v1, Table 1; H2-1E specification values.
# Depolarizing-parameter convention UNVERIFIED (process vs average infidelity).
# Use the requested numbers directly as total nonidentity Pauli probabilities.
H2_LIKE = IonParams(
    name='H2_LIKE',
    p_2q=1.28e-3,  # Table 1, two-qubit depolarizing error.
    p_1q=2.9e-5,  # Table 1, single-qubit depolarizing error.
    p_reset=4.0e-5,  # Table 1, initialization bit-flip.
    p_meas=1.0e-3,  # Table 1, measurement bit-flip.
    T2=None,  # Not supplied by this preset; no idle model in Stage 1.
    ms_duration=None,  # Not supplied; do not borrow another device's timing.
    rotation_duration=None,  # Not supplied.
    measurement_duration=None,  # Not supplied.
    reset_duration=None,  # Not supplied.
    shuttle_duration=None,  # Not supplied.
    split_duration=None,  # Not supplied.
    merge_duration=None,  # Not supplied.
    junction_duration=None,  # Not supplied.
    source='https://arxiv.org/html/2508.14227v1#S3.T1',
    caveat='Depolarizing-parameter convention unverified (process vs average infidelity). Timing unspecified.',
)

PRESETS = {p.name: p for p in (JONES_MURALI_2026, H2_LIKE)}


def fidelity_to_depolarizing(f_avg, n_qubits):
    """Average gate fidelity -> total probability of a nonidentity Pauli.

    For Hilbert dimension d, F_avg = 1 - d*p/(d+1). This is NOT the
    probability of replacement by the maximally mixed state.
    """
    if n_qubits not in (1, 2):
        raise ValueError('Only one- and two-qubit channels are supported')
    if not math.isfinite(f_avg) or not 0 <= f_avg <= 1:
        raise ValueError('Average fidelity must be in [0, 1]')
    d = 2 ** n_qubits
    p = (d + 1) / d * (1 - f_avg)
    if p > 1:
        raise ValueError('Fidelity is outside the range of this Pauli channel')
    return p

# Sensitivity values explicitly requested in the staged study prompt, in seconds.
T2_SWEEP = (0.1, 0.5, 2.2, 10.0)

# Stage 3 counterfactual readout-duration study, explicitly requested by the user.
# Error probabilities and all other nominal parameters remain fixed.
MEASUREMENT_DURATION_SWEEP = (0.0, 100e-6, 200e-6, 400e-6)

# Stage 3 analysis grid (not a device calibration). Two six-point neighborhoods
# share 0.00135; the union brackets H2-like, Jones–Murali/readout targets and
# the expected uniform MC pseudothreshold. Extend rather than extrapolate.
UNIFORM_CALIBRATION_P = (0.00085, 0.00095, 0.00105, 0.00115, 0.00125, 0.00135,
                         0.00145, 0.00155, 0.00165, 0.00175, 0.00185)

# Stage 4 placeholders: additive *probability* penalty per transport event,
# not a calibrated heating model. Zero preserves the prior operation channel.
DEFAULT_HEATING_PENALTY = dict.fromkeys(('shuttle', 'split', 'merge', 'junction'), 0.0)
# Synthetic sensitivity values requested to test monotonicity, not hardware data.
SYNTHETIC_HEATING_PENALTIES = (0.0, 5e-4, 1e-3)

# Design-choice experiment inputs explicitly requested in the experiment prompt.
# Analysis scan values, not additional hardware calibrations.
DESIGN_UNIFORM_P = (5e-4, 1e-3, 1.5e-3)
DESIGN_TWO_QUBIT_P = (5e-4, 1e-3, 2e-3)
# Stage A2 analysis-only calibration extension (scan inputs, not device values).
# Approved Stage 3 grid is reused; these points bracket lower-error targets.
DESIGN_CALIBRATION_EXTENSION_P = (0.00015, 0.00025, 0.00035, 0.00045,
                                  0.00055, 0.00065, 0.00075)
# Stage A3 analysis-only memory lengths (requested round counts, not hardware).
MEMORY_ROUNDS = (1, 2, 4, 8, 16)
