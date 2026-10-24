"""Stage A1 checkpoint only. A2 order Monte Carlo comparisons are not implemented here."""

import argparse
import ast
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import random

from flag_design import CHECK_ORDERS, FlagDesign, TWO_QUBIT_ONLY
from ion_noise import IdleIonNoise, IonNoise
from ion_params import (DESIGN_TWO_QUBIT_P, DESIGN_UNIFORM_P, H2_LIKE, JONES_MURALI_2026)
from steane_noise import UniformNoise, all_locations, validated_channel

ROOT = Path(__file__).resolve().parents[1]


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    temp.replace(path)


def source_fingerprint():
    names=('flag_design.py','steane_ft.py','steane_noise.py','ion_noise.py','ion_timing.py','ion_transport.py','ion_params.py')
    sources=[(ROOT/'src'/p).read_bytes() for p in names]
    # A2 appends analysis-only grid inputs; they do not affect any A1 experiment.
    # Keep the approved A1 identity for exactly its original parameter prefix.
    before,marker,after=sources[-1].partition(b'# Stage A2 analysis-only calibration extension')
    if marker:
        # Exclude only a literal analysis-grid assignment, never runtime overrides.
        try:
            statements=ast.parse(marker+after).body
            assignment=statements[0]
            value=ast.literal_eval(assignment.value)
            if (len(statements)==1 and isinstance(assignment,ast.Assign)
                    and len(assignment.targets)==1 and isinstance(assignment.targets[0],ast.Name)
                    and assignment.targets[0].id=='DESIGN_CALIBRATION_EXTENSION_P'
                    and isinstance(value,tuple) and all(type(p) is float and 0<p<1 for p in value)):
                sources[-1]=before
        except (SyntaxError,ValueError,TypeError,AttributeError,IndexError):
            pass
    return hashlib.sha256(b'\n'.join(sources)).hexdigest()


def frozen_artifacts():
    expected=json.loads((ROOT/'tests/fixtures/stages0_4_artifacts.json').read_text())
    actual={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in expected['sha256']}
    return {'frozen_commit':expected['frozen_commit'],'artifact_count':len(actual),
            'all_byte_identical':actual==expected['sha256'],
            'mismatches':[name for name in actual if actual[name]!=expected['sha256'][name]]}


def cases():
    rows=[(f'UniformNoise({p:g})', UniformNoise(p), 'timing-free') for p in DESIGN_UNIFORM_P]
    rows += [(f'TWO_QUBIT_ONLY({p:g})', TWO_QUBIT_ONLY(p), 'timing-free') for p in DESIGN_TWO_QUBIT_P]
    rows += [('H2_LIKE',IonNoise(H2_LIKE),'timing-free'),
             ('JONES_MURALI_2026',IdleIonNoise(JONES_MURALI_2026),'SEQ'),
             ('JONES_MURALI_2026',IdleIonNoise(JONES_MURALI_2026),'PAR_MR')]
    return rows


def geometry_signature(design):
    # Positive support, not its probability, determines an exhaustive fault test.
    data={'order':design.round_order,'locations':all_locations(design.circuits),
          'support':[[v for v,p in validated_channel(design,l).items() if p>0]
                     for l in all_locations(design.circuits)]}
    return hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()


def fault_matrix():
    output=ROOT/'results/stageA1_fault_checks.json'
    result={'stage':'A1','source_fingerprint':source_fingerprint(),'rows':[],'frozen_artifacts':frozen_artifacts()}
    if not result['frozen_artifacts']['all_byte_identical']:
        raise RuntimeError('Frozen artifacts changed before A1')
    reuse={}
    for order_name, order in CHECK_ORDERS.items():
        for preset, model, policy in cases():
            design=FlagDesign(model,order,timing_policy='SEQ' if policy=='timing-free' else policy)
            signature=geometry_signature(design)
            repeated=signature in reuse
            checks=reuse[signature] if repeated else design.check_faults()
            reuse[signature]=checks
            row={'order_name':order_name,'round_order':order,'preset':preset,'timing_policy':policy,
                 'geometry_signature':signature,'reused_identical_support_check':repeated,'checks':checks}
            result['rows'].append(row)
            result['unique_exhaustive_checks']=len(reuse)
            save(output,result)
            print(f'{order_name} / {preset} / {policy}: flagged {checks["flagged_fault_failures"]}/{checks["flagged_fault_runs"]}, '
                  f'bare-only {checks["bare_only_fault_failures"]}/{checks["bare_only_fault_runs"]}, '
                  f'incoming {checks["incoming_failures"]}/{checks["incoming_runs"]}, mismatches {checks["mismatches"]}',flush=True)
            if checks['single_fault_failures'] or checks['mismatches']:
                with (ROOT/'docs/STAGE_LOG.md').open('a') as fh:
                    fh.write('\n\n## Stage A1 — Fault-check failure; stopped\n\n')
                    fh.write(json.dumps(row,indent=2)+'\n')
                raise SystemExit('FT or simulator cross-check failure. Details saved; stop for review.')
    return result


def reference(policy, shots=5_000_000):
    if shots < 5_000_000:
        raise ValueError('Use at least 5M shots per reference arm')
    output=ROOT/f'results/stageA1_reference_{policy}.json'
    seed={'SEQ':73001,'PAR_MR':73002}[policy]
    config={'source_fingerprint':source_fingerprint(),'policy':policy,'shots':shots,'seed':seed,
            'parameters':asdict(JONES_MURALI_2026),'round_order':CHECK_ORDERS['ZZZXXX']}
    if output.exists():
        cached=json.loads(output.read_text())
        # JSON normalizes tuple inputs to lists.
        if cached['configuration']!=json.loads(json.dumps(config)):
            raise ValueError('Reference cache configuration changed; use a fresh output namespace')
        if cached.get('complete'):
            return cached
    design=FlagDesign(IdleIonNoise(JONES_MURALI_2026),round_order=CHECK_ORDERS['ZZZXXX'],timing_policy=policy)
    checks=design.check_faults()
    row={'configuration':config,'checks':checks,'complete':False}
    save(output,row)
    if checks['single_fault_failures'] or checks['mismatches']:
        raise SystemExit('Reference fault/cross-check failed; stop for review')
    row['counts']=design.count_faults()
    print(f'Baseline-order {policy}: P1={row["counts"]["P1"]}, P2={row["counts"]["P2"]}; {shots} MC shots',flush=True)
    row['monte_carlo']=design.monte_carlo(shots,random.Random(seed),return_details=True)
    row['monte_carlo']['seed']=seed
    row['idle_exposure']=design.exposure(row['monte_carlo'])
    row['complete']=True
    save(output,row)
    print(f'{policy}: MC={row["monte_carlo"]["p_L"]} +/- {row["monte_carlo"]["stderr"]}',flush=True)
    return row


def report():
    matrix=json.loads((ROOT/'results/stageA1_fault_checks.json').read_text())
    if matrix['source_fingerprint']!=source_fingerprint() or len(matrix['rows'])!=36:
        raise ValueError('Incomplete or stale A1 matrix')
    refs={policy:json.loads((ROOT/f'results/stageA1_reference_{policy}.json').read_text()) for policy in ('SEQ','PAR_MR')}
    if any(not r.get('complete') or r['configuration']['source_fingerprint']!=source_fingerprint() for r in refs.values()):
        raise ValueError('Incomplete or stale reference data')
    checks=[row['checks'] for row in matrix['rows']]+[r['checks'] for r in refs.values()]
    if any(c['single_fault_failures'] or c['mismatches'] for c in checks):
        raise ValueError('Fault or simulator check failed; stop for review')
    combined=math.hypot(*(r['monte_carlo']['stderr'] for r in refs.values()))
    precision=combined/min(r['monte_carlo']['p_L'] for r in refs.values())
    if precision>0.02:
        raise ValueError('Reference precision exceeds 2%; increase both shot counts')
    data={'stage':'A1','scope':'Order parameter and correctness only; A2/B not started.',
          'primary_estimate':'Monte Carlo','fault_matrix':matrix,'baseline_order_references':refs,
          'combined_reference_stderr':combined,'combined_stderr_fraction_of_smaller_p_L':precision,
          'reference_comparison_scope':'Timing-policy reference of the same baseline gadget, not an order-design comparison.',
          'timing_policy_contrast_only': {
              'delta_p_L_PAR_MR_minus_SEQ':refs['PAR_MR']['monte_carlo']['p_L']-refs['SEQ']['monte_carlo']['p_L'],
              'combined_stderr':combined,
              'sigma':(refs['PAR_MR']['monte_carlo']['p_L']-refs['SEQ']['monte_carlo']['p_L'])/combined,
              'delta_P2':refs['PAR_MR']['counts']['P2']-refs['SEQ']['counts']['P2']},
          'pair_attribution_vs_MC': {
              policy:{key:{'estimate_over_MC':r['counts'][key]/r['monte_carlo']['p_L'],
                           'estimate_minus_MC_in_stderr':(r['counts'][key]-r['monte_carlo']['p_L'])/r['monte_carlo']['stderr']}
                      for key in ('P2','P2_survival_corrected','P2_survival_mean_approx')}
              for policy,r in refs.items()},
          'frozen_artifacts':frozen_artifacts()}
    assert data['frozen_artifacts']['all_byte_identical']
    save(ROOT/'results/stageA1_validation.json',data)
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--part',choices=('checks','reference','report','all'),default='all')
    parser.add_argument('--policy',choices=('SEQ','PAR_MR'))
    parser.add_argument('--shots',type=int,default=5_000_000)
    args=parser.parse_args()
    if args.part in ('checks','all'): fault_matrix()
    if args.part in ('reference','all'):
        for policy in [args.policy] if args.policy else ('SEQ','PAR_MR'): reference(policy,args.shots)
    if args.part in ('report','all'): report()
