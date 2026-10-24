import json
import math
import numpy as np
import pytest
import stim

import steane_ft as ft
import steane_noise as old
from check_order import cases
from design_mc import BatchSampler
from flag_design import CHECK_ORDERS,FlagDesign,ec_round,ec_round_frame
from ion_noise import IdleIonNoise,_state_signature,_frame_signature
from ion_params import JONES_MURALI_2026
from memory_mc import (FlagMemoryRound,RoundStep,flag_endpoint_score,memory_trial,
                       merge_memory_counts,memory_metrics)
from memory_fit import decay_probability,fit_epsilon,compare_fits


def scalar_faults(sampler,sparse):
    result={}
    for li,(hits,codes) in sparse.items():
        loc=sampler.locations[li];variants=list(sampler.channels[li])
        for shot,code in zip(hits,codes):
            result.setdefault(shot,{}).setdefault(loc[:2],{})[loc[2]]=variants[code-1]
    return result


@pytest.mark.parametrize('order',['ZZZXXX','XXXZZZ','ZXZXZX'])
@pytest.mark.parametrize('policy',['SEQ','PAR_MR'])
def test_memory_multi_faults_carried_frame_scalar_and_both_basis_tableau(order,policy):
    d=FlagDesign(IdleIonNoise(JONES_MURALI_2026),CHECK_ORDERS[order],timing_policy=policy)
    a=FlagMemoryRound(d);rng=np.random.default_rng(18826);n=90;R=4
    faults=[]
    for _ in range(R):
        roundfaults={}
        for li,ch in enumerate(a.sampler.channels):
            codes=np.flatnonzero(np.array(list(ch.values()))>0)+1
            if len(codes):
                hits=np.flatnonzero(rng.random(n)<.012)
                roundfaults[li]=(hits,rng.choice(codes,len(hits)))
        faults.append(roundfaults)
    count,trace=memory_trial(a,n,R,rng,fault_rounds=faults,trace=True)
    per=[scalar_faults(a.sampler,f) for f in faults]
    for shot in range(n):
        x=z=0;sims={basis:stim.TableauSimulator() for basis in ('0','+')}
        for basis,sim in sims.items():ft.prepare_logical(sim,basis)
        for r in range(R):
            incoming=(x,z)
            log={};x,z=ec_round_frame(x,z,per[r].get(shot,{}),circuits=d.circuits,round_order=d.round_order,log=log)
            out=trace['rounds'][r]
            assert (x,z)==(out['state'][0][shot],out['state'][1][shot])
            assert bool(log['trigger'])==out['branch_taken'][shot]
            if log['trigger']:
                t=out['diagnostics']['trace'];j=np.searchsorted(t['shot_indices'],shot)
                expected=[]
                for category in ('reset','measurement','2q','idle','carry'):
                    subset={}
                    if category!='carry':
                        for key,hits in per[r].get(shot,{}).items():
                            for idx,payload in hits.items():
                                kind=old.loc_type((*key,idx,d.circuits[key][idx]))
                                kind='measurement' if kind=='meas' else kind
                                if kind==category:subset.setdefault(key,{})[idx]=payload
                    xx,zz=incoming if category=='carry' else (0,0)
                    for pos in range(log['trigger_index']+1):
                        xx,zz,bits=old.run_ops_frame(xx,zz,d.circuits[('F',pos)],subset.get(('F',pos),{}))
                    expected.append(any(bits[tag] and log['trigger_bits'][tag] for tag in ('synd','flag')))
                assert expected==t['contributors'][:,j].tolist()
            for basis,sim in sims.items():
                assert ec_round(sim,faults=per[r].get(shot,{}),circuits=d.circuits,round_order=d.round_order)==log
                assert _state_signature(sim,basis)==_frame_signature(x,z,basis)
        assert old.logical_failure(x,z)==(trace['fail_X'][shot],trace['fail_Z'][shot])
    assert count['round_totals']['triggers']==sum(c['triggers'] for c in count['per_round'])


@pytest.mark.parametrize('order',['ZZZXXX','XXXZZZ','ZXZXZX'])
def test_memory_R1_exactly_matches_frozen_A2_explicit_shots(order):
    for _,noise,policy in cases():
        d=FlagDesign(noise,CHECK_ORDERS[order],timing_policy='SEQ' if policy=='timing-free' else policy)
        a=FlagMemoryRound(d);rng=np.random.default_rng(4426)
        faults=a.sampler.sample_faults(2000,rng)
        frozen=BatchSampler(d).simulate(2000,faults)
        counts=memory_trial(a,2000,1,rng,fault_rounds=[faults])
        assert counts['fails']==frozen['fails']
        assert counts['logical_failure_counts']==frozen['logical_failure_counts']
        assert counts['round_totals']['triggers']==frozen['triggers']
        for c in ('reset','measurement','2q','idle','mixed'):
            assert counts['round_totals']['trigger_cause_counts'][c]==frozen['trigger_cause_counts'][c]
        assert counts['round_totals']['trigger_cause_counts']['carry']==0


def test_no_intermediate_ideal_decode_and_generic_adapter_abort_contract():
    # A toy adapter leaves X0 in round one, then adds X1. An intermediate ideal
    # correction would erase X0 and incorrectly make the final shot successful.
    class Toy:
        def __init__(self):self.r=0
        def initial_state(self,n):return np.zeros(n,np.uint16),np.zeros(n,np.uint16)
        def step(self,state,rng,*,faults=None):
            x,z=state;out=x^(1<<self.r);self.r+=1
            return RoundStep((out,z.copy()),np.zeros(len(x),bool),np.zeros(len(x),bool),{})
    count=memory_trial(Toy(),5,2,np.random.default_rng(0))
    assert count['fails_X']==5 and count['fails']==5
    class LogicalFlips(Toy):
        def step(self,state,rng,*,faults=None):
            x,z=state
            return RoundStep((x^7,z.copy()),np.zeros(len(x),bool),np.zeros(len(x),bool),{})
    assert memory_trial(LogicalFlips(),5,1,np.random.default_rng(0))['fails']==5
    assert memory_trial(LogicalFlips(),5,2,np.random.default_rng(0))['fails']==0
    class Aborts:
        def initial_state(self,n):return np.zeros(n,np.uint16),np.zeros(n,np.uint16)
        def step(self,state,rng,*,faults=None):
            n=len(state[0]);return RoundStep(state,np.zeros(n,bool),np.ones(n,bool),{})
    assert memory_trial(Aborts(),3,2,np.random.default_rng(0))['aborted']==3


def test_noiseless_round_cleans_incoming_weight_one_without_resetting_logical():
    a=FlagMemoryRound(FlagDesign(old.UniformNoise(0)))
    for p in 'XYZ':
        x=np.array([1<<q if p in 'XY' else 0 for q in range(7)],np.uint16)
        z=np.array([1<<q if p in 'ZY' else 0 for q in range(7)],np.uint16)
        step=a.step((x,z),np.random.default_rng(1),faults={})
        assert not any(np.any(c) for c in flag_endpoint_score(step.state))
        assert step.branch_taken.all()
        assert step.diagnostics['trigger_cause_counts']['carry']==7
    # Logical X on a weight-3 representative has zero syndrome and persists.
    logical=7
    step=a.step((np.array([logical],np.uint16),np.zeros(1,np.uint16)),np.random.default_rng(1),faults={})
    assert flag_endpoint_score(step.state)[0][0]


def test_fit_known_decay_excludes_R1_and_flags_boundary_misspecification():
    epsilon=.002;n=10_000_000
    points=[{'rounds':r,'shots':n,'fails':round(n*float(decay_probability(epsilon,r)))} for r in (1,2,4,8,16)]
    points[0]['fails']=n//5
    fit=fit_epsilon(points)
    assert fit['epsilon']==pytest.approx(epsilon,rel=1e-5)
    assert fit['relative_stderr']<.03 and fit['goodness_of_fit_p_value']>.9
    biased=[dict(p,fails=max(1,p['fails']-10000)) for p in points]
    assert fit_epsilon(biased)['goodness_of_fit_p_value']<1e-8
    c=compare_fits(fit,dict(fit,epsilon=fit['epsilon']+.0001))
    assert c['delta_epsilon']<0 and c['significant_at_3sigma']


def test_memory_sampling_merge_and_resume_state():
    a=FlagMemoryRound(FlagDesign(old.UniformNoise(.01)));rng=np.random.default_rng(7726)
    first=memory_trial(a,500,4,rng);state=json.loads(json.dumps(rng.bit_generator.state))
    second=memory_trial(a,500,4,rng)
    rng2=np.random.default_rng();rng2.bit_generator.state=state
    assert second==memory_trial(a,500,4,rng2)
    merged=merge_memory_counts(first,second);metrics=memory_metrics(merged)
    assert metrics['shots']==1000
    assert sum(metrics['round_totals']['trigger_cause_counts'].values())==metrics['round_totals']['triggers']
    assert sum(metrics['logical_failure_counts'].values())==metrics['fails']


@pytest.mark.parametrize('order',['ZZZXXX','XXXZZZ','ZXZXZX'])
def test_every_supported_single_fault_early_middle_last_of_R16(order):
    from memory_mc import memory_single_fault_checks
    from ion_noise import IonNoise
    from ion_params import H2_LIKE
    for noise,policy in [(old.UniformNoise(.001),'SEQ'),(IonNoise(H2_LIKE),'SEQ'),
                         (IdleIonNoise(JONES_MURALI_2026),'SEQ'),(IdleIonNoise(JONES_MURALI_2026),'PAR_MR')]:
        a=FlagMemoryRound(FlagDesign(noise,CHECK_ORDERS[order],timing_policy=policy))
        checks=memory_single_fault_checks(a)
        assert checks['payload_runs']>0 and checks['failures_union']==0
        assert checks['failures_X']==checks['failures_Z']==0


def test_memory_cache_resumption_and_configuration_guard(tmp_path):
    from memory_experiment import cache_point
    a=FlagMemoryRound(FlagDesign(old.UniformNoise(.01)))
    path=tmp_path/'memory.json'
    first=cache_point(a,4,500,8926,path)
    final=cache_point(a,4,1000,8926,path)
    assert final['chunks']==2 and final['counts']['shots']==1000
    rng=np.random.default_rng(8926)
    expected=merge_memory_counts(memory_trial(a,500,4,rng),memory_trial(a,500,4,rng))
    assert final['counts']==expected
    assert cache_point(a,4,1000,8926,path)==final
    with pytest.raises(ValueError,match='incompatible'):cache_point(a,8,1000,8926,path)
    with pytest.raises(ValueError,match='incompatible'):cache_point(a,4,1000,8927,path)


def test_exact_pair_mechanism_conserves_excess_and_matches_replays():
    from pathlib import Path
    import hashlib
    root=Path(__file__).parents[1]
    r=json.loads((root/'results/stageA3_pair_mechanism.json').read_text())
    assert r['source_sha256']==hashlib.sha256((root/'src/pair_mechanism.py').read_bytes()).hexdigest()
    parts=r['matched_payload_partition']
    assert parts['new_only']['P2_weight']-parts['removed']['P2_weight']==pytest.approx(r['delta_P2'])
    assert sum(c['delta_P2'] for c in r['by_location_types'])==pytest.approx(r['delta_P2'])
    classes=r['classes_positive_delta']+r['classes_negative_delta']
    assert sum(c['delta_P2'] for c in classes)==pytest.approx(r['delta_P2'])
    assert r['delta_P2_survival_corrected']==pytest.approx(r['delta_P2']*(1-r['p'])**94)
    for example in r['top_examples']:
        p=example['fault_propagation']
        assert not any(p['baseline_replay']['logical_failure'])
        assert any(p['alternating_replay']['logical_failure'])


def test_combined_basis_metric_propagates_shared_shot_covariance():
    from memory_fit import combine_basis_fits
    n=10_000_000;points=[]
    for r in (1,2,4,8,16):
        px=float(decay_probability(.001,r));pz=float(decay_probability(.002,r))
        kx=round(n*px);kz=round(n*pz);both=round(n*.15*min(px,pz))
        points.append({'rounds':r,'shots':n,'fails_X':kx,'fails_Z':kz,
                       'logical_failure_counts':{'both_XZ':both}})
    fx=fit_epsilon(points,metric='fails_X');fz=fit_epsilon(points,metric='fails_Z')
    total=combine_basis_fits(points,fx,fz)
    assert total['epsilon']==pytest.approx(1-(1-fx['epsilon'])*(1-fz['epsilon']))
    assert total['covariance_epsilon_X_Z']>0
    assert total['stderr']>total['stderr_if_covariance_ignored_secondary']
    assert total['relative_stderr']<.03 and not total['union_is_fitted']
    # Check the MLE influence against a direct perturbation of an endpoint count.
    changed=[dict(p) for p in points];changed[1]['fails_X']+=10
    changed_fit=fit_epsilon(changed,metric='fails_X')
    influence=fx['points'][0]['epsilon_influence_per_failure_count']
    assert changed_fit['epsilon']-fx['epsilon']==pytest.approx(10*influence,rel=1e-4)


def test_R2_deviation_is_checked_against_held_out_tail_trend():
    from memory_fit import R2_trend_diagnostic
    n=1_000_000
    points=[{'rounds':r,'shots':n,'fails_X':round(n*float(decay_probability(.001,r)))} for r in (1,2,4,8,16)]
    points[1]['fails_X']+=300
    tail=fit_epsilon(points,metric='fails_X',minimum_rounds=4)
    diagnostic=R2_trend_diagnostic(points,tail,'fails_X')
    assert diagnostic['outside_R2_one_sigma_MC_errorbar']
    assert diagnostic['sigma_combined']>3 and tail['fit_rounds']==[4,8,16]
