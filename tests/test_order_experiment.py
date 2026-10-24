import json
from pathlib import Path
import numpy as np
import pytest

from flag_design import FlagDesign
from order_experiment import (_scale_pair_counts, frozen_a1, sample_cached,
                               simulation_fingerprint)
from steane_noise import UniformNoise

ROOT=Path(__file__).parents[1]


def test_pair_reweighting_matches_approved_uniform_coefficients():
    c=json.loads((ROOT/'results/stage0_channels.json').read_text())
    for p in (0.0005,0.001,0.0015):
        # Original Stage 0 counts were evaluated at p=.0005.
        scaled=_scale_pair_counts(c,p,0.0005,96)
        assert scaled['P2']==pytest.approx(c['A']['any']*p*p)
        assert scaled['P2_survival_corrected']==pytest.approx(scaled['P2']*(1-p)**94)
        assert sum(scaled['P2_by_type_pair'].values())==pytest.approx(scaled['P2'])
        assert sum(scaled['P2_survival_by_type_pair'].values())==pytest.approx(scaled['P2_survival_corrected'])
        assert scaled['P2_by_logical_type']['any']==pytest.approx(scaled['P2'])


def test_approved_a1_artifacts_unchanged():
    assert len(frozen_a1())==4 and all(frozen_a1().values())


def test_mc_chunk_cache_resume_and_configuration_guard(tmp_path):
    d=FlagDesign(UniformNoise(.01));path=tmp_path/'mc.json'
    first=sample_cached(path,d,2000,413)
    resumed=sample_cached(path,d,4000,413)
    assert resumed['shots']==4000
    data=json.loads(path.read_text());assert data['chunks']==2
    assert data['configuration']['simulation_fingerprint']==simulation_fingerprint()
    assert sample_cached(path,d,4000,413)==resumed
    # Replay precisely the same two chunks from an empty cache.
    other=tmp_path/'other.json'
    assert sample_cached(other,d,2000,413)==first
    assert sample_cached(other,d,4000,413)==resumed
    with pytest.raises(ValueError,match='incompatible'):sample_cached(path,d,5000,414)
