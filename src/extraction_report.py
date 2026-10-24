"""Assemble the B2 resource/checkpoint report; no performance calculations."""
import json
import math
from collections import Counter
from pathlib import Path
from check_order import save
from extraction_experiment import designs,fingerprint,frozen_checks
from prep_experiment import bernoulli

ROOT=Path(__file__).resolve().parents[1]


def assemble():
    checked=json.loads((ROOT/'results/stageB2_correctness.json').read_text())
    smoke=json.loads((ROOT/'results/stageB2_memory_smoke.json').read_text())
    assert checked['complete'] and smoke['complete']
    assert checked['simulation_fingerprint']==smoke['simulation_fingerprint']==fingerprint()
    b1=json.loads((ROOT/'results/stageB1_preparation.json').read_text())
    output={'stage':'B2','scope':'Correctness, resources and harness validation only. No B3 performance estimate or comparison.',
        'simulation_fingerprint':fingerprint(),'harness_fingerprint':smoke['harness_fingerprint'],'approved_artifacts':frozen_checks(),
        'max_attempts':4,'ion_count':15,'rows':{},
        'deferred_to_B3':['epsilon_tot and its abort share','epsilon_X and epsilon_Z performance fits','Jones–Murali SEQ max_attempts in {2,4,infinity} performance sensitivity','design comparisons and calibration'],
        'performance_metrics':{'epsilon_tot':None,'epsilon_X':None,'epsilon_Z':None,'abort_share_of_epsilon_tot':None}}
    for key,preset,policy,d in designs():
        row=checked['rows'][key];ci=key.split('_')[0]
        qs=[b1['rows'][ci+('_plus' if kind=='X' else '_zero')]['exact_finite_rate']['rejection'] for kind in d.order]
        weights={};reach=1.
        for hi,q in enumerate(qs):
            for a in range(4):weights[('P',hi,a)]=reach*q**a
            reach*=1-q**4;weights[('C',hi,0)]=reach
        expected_counts=Counter();expected_duration=0.;exposure={str(q):{'seconds':0.,'summed_Z_probability':0.} for q in range(7)}
        for segment,w in weights.items():
            resource=d.resources([segment])
            for typ,n in resource['location_counts'].items():expected_counts[typ]+=w*n
            expected_duration+=w*(resource['duration_seconds'] or 0.)
            if d.timed:
                for q in exposure:
                    for field in exposure[q]:exposure[q][field]+=w*resource['idle_exposure_per_data_qubit'][q][field]
        c=row['checks']['counts'];nominal_n=sum(c[k] for k in ('noiseless_runs','first_attempt_or_extraction_fault_runs','retry_only_fault_runs','incoming_single_error_runs'))
        r0,r1=qs[0]**4,qs[1]**4
        output['rows'][key]={'preset':preset,'policy':policy,'half_order':d.half_order,'half_sequence':d.order,
            'no_rejection':row['no_rejection_resources'],'maximum_executed_if_fourth_attempt_accepts_in_each_half':row['all_reserved_locations'],
            'model_expected_resources_including_rejections_and_aborts':{'location_counts':dict(expected_counts),
                'duration_seconds':expected_duration if d.timed else None,'idle_exposure_per_data_qubit':exposure if d.timed else None,
                'method':'Per-attempt verification rejection from B1 exact finite-rate convolution; actual execution reach weights q^a and previous-half acceptance. This is a resource expectation, not p_L.'},
            'model_expected_round_abort_probability':r0+r1-r0*r1,
            'observed_unconditioned_correctness_aborts':{'aborts':0,'round_runs':nominal_n,'rate':0.,'scope':'Deterministic fault-injection/noiseless/incoming-error trials; not a stochastic physical abort-rate estimate.'},
            'observed_conditioned_retry_aborts':{'aborts':c['expected_conditioned_aborts'],'round_runs':c['conditioned_retry_fault_runs'],'scope':'Test-only branch-conditioned retry paths; not physical sampling.'},
            'physical_exhaustion_tests':{'aborts':c['physical_abort_policy_runs'],'round_runs':c['physical_abort_policy_runs'],'counts_as_failure':True}}
    n=aborts=0
    for row in smoke['rows'].values():
        for entry in row['smoke']:
            n+=entry['executed_rounds'];aborts+=entry['observed_aborts']['count']
    output['observed_smoke_abort_rate']=dict(bernoulli(aborts,n),scope='Pooled actual executed Jones–Murali rounds across both half-orders/policies; harness sanity only, not a performance comparison.')
    save(ROOT/'results/stageB2_resources.json',output)
    return output


if __name__=='__main__':assemble()
