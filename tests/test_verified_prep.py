import itertools
import json
import math
from pathlib import Path

import numpy as np
import pytest
import stim

import steane_ft as ft
import steane_noise as scalar
from flag_design import layer_timing
from ion_noise import IonNoise,IdleIonNoise
from ion_params import H2_LIKE,JONES_MURALI_2026
from ion_timing import operation_duration
from prep_experiment import derivation,frozen_checks,sample_cache,DEFAULT_MAX_ATTEMPTS
from verified_prep import (VerifiedPreparation,BLOCK,VERIFIER,PAPER_ENCODER,VERIFY_SUPPORT,
                           CODE_GROUP,REDUCED_WEIGHT,STATE_PHASE_WEIGHT,preparation_ops,
                           run_tableau,tableau_signature,frame_signature,weights,popcount)

ROOT=Path(__file__).parents[1]


def test_published_encoder_row_span_and_weight_three_verification_derivation():
    d=derivation()
    assert d['same_stabilizer_row_span_as_frozen_code']
    assert d['generated_X_stabilizer_supports']==[[0,1,4,5],[0,2,4,6],[3,4,5,6]]
    assert sum(v['payload_count'] for v in d['unverified_dangerous_cosets'].values())==30
    assert len(d['unverified_dangerous_cosets'])==2
    assert all(v['anticommutes_with_verifier'] for v in d['unverified_dangerous_cosets'].values())
    assert VERIFY_SUPPORT==(0,5,6)
    a,b=preparation_ops('0'),preparation_ops('+')
    assert len(a)==len(b)==20 and sum(op[0]=='CX' for op in a)==11
    for x,y in zip(a,b):
        if x[0]=='CX':assert y==('CX',x[2],x[1])
        elif x[0] in ('R','RX'):assert y[0]==('RX' if x[0]=='R' else 'R') and y[1]==x[1]
        else:assert x==('M',VERIFIER,'verify') and y==('MX',VERIFIER,'verify')


@pytest.mark.parametrize('basis',['0','+'])
@pytest.mark.parametrize('policy',['timing-free','SEQ','PAR_MR'])
def test_every_universal_single_fault_rejected_or_safe_and_tableau_agrees(basis,policy):
    noise=scalar.UniformNoise(.001) if policy=='timing-free' else IdleIonNoise(JONES_MURALI_2026)
    g=VerifiedPreparation(basis,noise,timing_policy='SEQ' if policy=='timing-free' else policy)
    assert g.noiseless_check()['verification_outcome']==0
    checks=g.check_faults(universal=True)
    assert checks['prep_fault_runs']>=190
    assert checks['accepted_dangerous_failures']==checks['frame_tableau_mismatches']==0
    assert checks['accepted_maxima']=={'other_raw':3,'other_code_reduced':3,
                                       'other_prepared_state_reduced':1,'dangerous_code_reduced':1}


def test_dangerous_logical_operator_is_not_removed_by_prepared_state_degeneracy():
    # X0 X1 X2 is logical X: harmful on |0_L>, stabilizing on |+_L>.
    x=7<<7
    assert weights(x,0,'0')['dangerous_code_reduced']==3
    assert weights(x,0,'+')['other_prepared_state_reduced']==0
    assert weights(x,0,'+')['other_code_reduced']==3
    assert weights(0,x,'+')['dangerous_code_reduced']==3
    assert weights(0,x,'0')['other_prepared_state_reduced']==0


@pytest.mark.parametrize('basis',['0','+'])
@pytest.mark.parametrize('policy',['SEQ','PAR_MR'])
def test_preparation_live_windows_conserve_time_and_include_verifier_readout_wait(basis,policy):
    params=JONES_MURALI_2026;g=VerifiedPreparation(basis,IdleIonNoise(params),timing_policy=policy)
    starts,ends,total,_=layer_timing(g.base_ops,params,policy)
    assert [op for op in g.ops if op[0]!='IDLE']==g.base_ops
    expected=(5*params.reset_duration+3*(params.reset_duration+params.rotation_duration)+11*params.ms_duration+params.measurement_duration
              if basis=='0' else 3*params.reset_duration+5*(params.reset_duration+params.rotation_duration)+11*params.ms_duration+params.rotation_duration+params.measurement_duration)
    if policy=='PAR_MR':expected=params.reset_duration+params.rotation_duration+11*params.ms_duration+params.measurement_duration+(params.rotation_duration if basis=='+' else 0)
    assert total==pytest.approx(expected) and g.duration==total
    for q in (*BLOCK,VERIFIER):
        reset=next(i for i,o in enumerate(g.base_ops) if o[0] in ('R','RX') and o[1]==q)
        busy=sum(operation_duration(o,params) for o in g.base_ops if o[0]=='CX' and q in o[1:])
        end=starts[19] if q==VERIFIER else total
        actual=sum(o[2] for o in g.ops if o[0]=='IDLE' and o[1]==q)
        assert actual==pytest.approx(end-ends[reset]-busy,abs=1e-15)
        if q in BLOCK:
            readout_wait=sum(o[2] for o,m in zip(g.ops,g.metadata) if o[0]=='IDLE' and o[1]==q and m['part']=='verification_readout')
            assert readout_wait==pytest.approx(operation_duration(g.base_ops[19],params))
    assert not any(o[0]=='IDLE' and o[1]<7 for o in g.ops)


@pytest.mark.parametrize('basis',['0','+'])
@pytest.mark.parametrize('policy',['timing-free','SEQ','PAR_MR'])
def test_full_frame_linear_influences_and_random_codeword_syndromes(basis,policy):
    noise=scalar.UniformNoise(.05) if policy=='timing-free' else IdleIonNoise(JONES_MURALI_2026)
    g=VerifiedPreparation(basis,noise,timing_policy='SEQ' if policy=='timing-free' else policy)
    rng=np.random.default_rng(55126);raw_words=set()
    for _ in range(80):
        faults={};x=z=verify=0
        for idx,(channel,influences) in enumerate(zip(g.channels,g.influences)):
            positive=[j for j,p in enumerate(channel.values()) if p>0]
            if positive and rng.random()<.07:
                j=int(rng.choice(positive));v=influences[j];faults[idx]=v['payload']
                x ^=v['fx'];z ^=v['fz'];verify ^=v['verify']
        assert (x,z,verify)==g.run_frame(faults)
        sim,bits=run_tableau(g.ops,faults)
        assert bits['verify']==verify
        assert tableau_signature(sim,basis)==frame_signature(x,z,basis)
        # Ideal terminal diagnostic readout only; not part of the B1 gadget.
        # It is the readout basis intended for later syndrome extraction.
        raw=[]
        for q in BLOCK:
            if basis=='0':sim.h(q)
            raw.append(int(sim.measure(q)))
        raw_words.add(tuple(raw))
        measured=[sum(raw[q] for q in support)%2 for support in ft.SUPPORTS]
        frame=((z if basis=='0' else x)>>7)&127
        assert measured==[popcount(frame&s)%2 for s in [sum(1<<q for q in support) for support in ft.SUPPORTS]]
    assert len(raw_words)>1 # raw words are deliberately not compared shot-by-shot


def test_exact_convolution_two_fault_survival_and_unequal_sampling():
    base=VerifiedPreparation('0',scalar.UniformNoise(.001));selected=None
    for a,b in itertools.combinations(range(8,len(base.ops)),2):
        for va in base.influences[a]:
            for vb in base.influences[b]:
                code=va['code']^vb['code']
                if code<128 and REDUCED_WEIGHT[code&127]>1:
                    selected=(a,va['payload'],b,vb['payload']);break
            if selected:break
        if selected:break
    assert selected is not None
    a,va,b,vb=selected
    class Sparse:
        def channel(self,loc):
            idx=loc[2]
            if idx==a:return {va:.1}
            if idx==b:return {vb:.2}
            if idx==7:return {'Z':.3} # verifier reset Z is inert for this objective
            return {}
    g=VerifiedPreparation('0',Sparse());exact=g.exact_distribution();pairs=g.count_pairs()
    assert exact['accepted_dangerous_probability']==pytest.approx(.02)
    assert pairs['P1']==0 and pairs['P2']==pytest.approx(.02)
    assert pairs['P2_survival_corrected']==pytest.approx(.014)
    counts=g.sample_attempts(100000,np.random.default_rng(98126))
    assert counts['accepted_dangerous']/counts['shots']==pytest.approx(.02,abs=.003)


def test_output_arrays_do_not_change_sampling_and_match_dangerous_component():
    g=VerifiedPreparation('+',IdleIonNoise(JONES_MURALI_2026),timing_policy='PAR_MR')
    expected=g.sample_attempts(10000,np.random.default_rng(87326))
    counts,arrays=g.sample_attempts(10000,np.random.default_rng(87326),return_arrays=True)
    assert counts==expected
    assert np.array_equal(arrays['codes']&127,(arrays['fz']>>7)&127)
    assert sum(counts['verification_rejection_cause_counts'].values())==counts['rejected']
    assert sum(counts['verification_fractional_type_counts'].values())==pytest.approx(counts['rejected'])


def test_verification_probability_has_independent_parity_product_benchmark():
    for basis in ('0','+'):
        for policy in ('SEQ','PAR_MR'):
            g=VerifiedPreparation(basis,IdleIonNoise(JONES_MURALI_2026),timing_policy=policy)
            flip=[sum(p for p,v in zip(ch.values(),table) if v['verify']) for ch,table in zip(g.channels,g.influences)]
            parity_accept=(1+math.prod(1-2*p for p in flip))/2
            assert g.exact_distribution()['acceptance']==pytest.approx(parity_accept,abs=2e-14)


def test_retry_attempt_count_exhaustion_and_fresh_resets():
    class AlwaysReject:
        def channel(self,loc):return {'FLIP':1.} if loc[3][0]=='M' else {}
    g=VerifiedPreparation('0',AlwaysReject())
    assert DEFAULT_MAX_ATTEMPTS==4
    result=g.sample_retries(20,np.random.default_rng(4),max_attempts=DEFAULT_MAX_ATTEMPTS)
    assert result['max_attempts']==4 and result['accepted_on_attempt']==[0,0,0,0]
    assert result['all_attempts_rejected']==20 and result['attempt_executions']==80
    with pytest.raises(ValueError):g.sample_retries(20,np.random.default_rng(4),max_attempts=0)
    # Fresh preparation resets all block/verifier errors, preserving spectators.
    clean=VerifiedPreparation('0',scalar.UniformNoise(0))
    x,z,verify=clean.run_frame(fx=(1<<15)-1,fz=(1<<15)-1)
    assert (x,z,verify)==(127,127,0)


def test_MC_cache_resumes_exactly_and_guards_seed(tmp_path):
    g=VerifiedPreparation('0',scalar.UniformNoise(.01));p=tmp_path/'mc.json'
    first=sample_cache(g,1000,92126,p);final=sample_cache(g,2000,92126,p)
    other=tmp_path/'other.json'
    assert sample_cache(g,1000,92126,other)==first
    assert sample_cache(g,2000,92126,other)==final
    with pytest.raises(ValueError,match='Stale'):sample_cache(g,3000,92127,p)


def test_timing_free_H2_and_transport_scope_and_approved_artifacts():
    assert frozen_checks()['artifact_count']==147
    assert frozen_checks()['all_byte_identical']
    with pytest.raises(ValueError,match='Untimed'):VerifiedPreparation('0',IonNoise(H2_LIKE),timing_policy='PAR_MR')
    moved=IdleIonNoise(JONES_MURALI_2026,transport_events=[('F',0,0,0,'shuttle')])
    with pytest.raises(ValueError,match='transport'):VerifiedPreparation('0',moved)
    assert all(op[0]!='IDLE' for op in VerifiedPreparation('+',IonNoise(H2_LIKE)).ops)
