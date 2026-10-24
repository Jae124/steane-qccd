"""Exact matched-physical-fault A3 mechanism study for Uniform(1e-3).

Compare identical (stage, kind, support, original op index, payload) fault
pairs in the two executed orders. Full positive/negative classes retain both
position maps and both trigger paths, so reordering is not mistaken for a
change of physical fault identity. No probabilities or circuits are adjusted.
"""
from collections import defaultdict
import itertools
import hashlib
import json
from pathlib import Path
import numpy as np

from check_order import save
from design_mc import BatchSampler
from flag_design import CHECK_ORDERS,FlagDesign,ec_round_frame
import steane_noise as scalar

ROOT=Path(__file__).resolve().parents[1]


def identity(design,loc):
    stage,pos,j,op=loc;kind,k=design.round_order[pos]
    return stage,kind,k,j


def _pauli(x,z,n=9):
    return ''.join('Y' if (x>>q)&1 and (z>>q)&1 else 'X' if (x>>q)&1 else 'Z' if (z>>q)&1 else 'I' for q in range(n))


def detailed_replay(design,faults):
    log={};expected=ec_round_frame(0,0,faults,circuits=design.circuits,round_order=design.round_order,log=log)
    x=z=0;events=[];prefix=[]
    keys=([('F',i) for i in range(6)] if log['trigger_index'] is None else
          [('F',i) for i in range(log['trigger_index']+1)]+[('B',i) for i in range(6)])
    for stage,pos in keys:
        kind,k=design.round_order[pos];prefix.append(f'{stage}{pos}:{kind}{k}')
        for j,op in enumerate(design.circuits[(stage,pos)]):
            payload=faults.get((stage,pos),{}).get(j);before=(x,z)
            x,z,bits=scalar.run_ops_frame(x,z,[op],{} if payload is None else {0:payload})
            if payload is not None or bits or before!=(x,z):
                events.append({'circuit':f'{stage}{pos}:{kind}{k}','op_index':j,'op':op,
                    'payload':payload,'frame_before':_pauli(*before),'frame_after':_pauli(x,z),'measurement':bits})
    before_correction=_pauli(x,z)
    for q in log.get('x_fix',[]):x ^=1<<q
    for q in log.get('z_fix',[]):z ^=1<<q
    assert (x,z)==expected
    return {'executed_checks':prefix,'log':log,'before_correction':before_correction,
            'residual_frame':_pauli(x,z),'logical_failure':list(scalar.logical_failure(x,z)),
            'events':events}


def run_mechanism():
    # Explicitly use the requested existing uniform scan value.
    from ion_params import DESIGN_UNIFORM_P
    p=DESIGN_UNIFORM_P[1]
    designs=[FlagDesign(scalar.UniformNoise(p),CHECK_ORDERS[name]) for name in ('ZZZXXX','ZXZXZX')]
    samplers=[BatchSampler(d) for d in designs];base,alt=samplers
    alternative={identity(designs[1],loc):i for i,loc in enumerate(alt.locations)}
    mapping=np.array([alternative[identity(designs[0],loc)] for loc in base.locations])
    for i,j in enumerate(mapping):
        assert base.locations[i][3]==alt.locations[j][3] and base.channels[i]==alt.channels[j]
    aa=[];bb=[];va=[];vb=[];weights=[]
    for a,b in itertools.combinations(range(len(base.locations)),2):
        na,nb=len(base.channels[a]),len(base.channels[b])
        count=na*nb
        aa.extend([a]*count);bb.extend([b]*count)
        va.extend(np.repeat(np.arange(1,na+1),nb));vb.extend(np.tile(np.arange(1,nb+1),na))
        weights.extend([p*p/(na*nb)]*count)
    aa=np.array(aa);bb=np.array(bb);va=np.array(va,np.uint8);vb=np.array(vb,np.uint8);w=np.array(weights)
    n=len(aa);traces=[];results=[]
    for sampler,locmap in [(base,np.arange(len(base.locations))),(alt,mapping)]:
        faults={}
        for li in range(len(base.locations)):
            chosen=np.flatnonzero((aa==li)|(bb==li))
            faults[int(locmap[li])]=(chosen,np.where(aa[chosen]==li,va[chosen],vb[chosen]))
        result=sampler.simulate(n,faults,traces=True);traces.append(result.pop('traces'));results.append(result)
    failed=[t['x_fail']|t['z_fail'] for t in traces]
    survival=(1-p)**(len(base.locations)-2)
    partition={label:{'payload_pairs':int(mask.sum()),'P2_weight':float(w[mask].sum()),
                      'P2_survival_weight':float(w[mask].sum()*survival)}
               for label,mask in [('new_only',failed[1]&~failed[0]),('removed',failed[0]&~failed[1]),
                                   ('common',failed[0]&failed[1])]}
    def path(i,shot):
        t=traces[i];pos=int(t['stop'][shot])
        if pos<0:return 'no_trigger'
        kind,k=designs[i].round_order[pos]
        return f'F{pos}:{kind}{k}:flag{int(t["flag"][shot])}:synd{int(t["trigger_synd"][shot])}->bare'
    def descriptor(li):
        l=base.locations[li];ident=identity(designs[0],l)
        return f'{ident[0]}:{ident[1]}{ident[2]}'
    groups={};types=defaultdict(lambda:np.zeros(2));positions=defaultdict(lambda:np.zeros(2))
    for shot in np.flatnonzero(failed[0]|failed[1]):
        a,b=int(aa[shot]),int(bb[shot]);la,lb=base.locations[a],base.locations[b]
        typepair='*'.join(sorted([scalar.loc_type(la),scalar.loc_type(lb)]))
        physical=(descriptor(a),descriptor(b));op_pair=(la[3][0],lb[3][0])
        bpos=(f'{la[0]}{la[1]}',f'{lb[0]}{lb[1]}')
        al,bl=alt.locations[mapping[a]],alt.locations[mapping[b]]
        apos=(f'{al[0]}{al[1]}',f'{bl[0]}{bl[1]}')
        paths=(path(0,shot),path(1,shot))
        key=(typepair,physical,op_pair,bpos,apos,paths)
        if key not in groups:
            groups[key]={'location_types':typepair,'physical_checks':physical,'operation_types':op_pair,
                         'baseline_positions':bpos,'alternating_positions':apos,
                         'baseline_trigger_path':paths[0],'alternating_trigger_path':paths[1],
                         'baseline_P2':0.,'alternating_P2':0.,'new_only_P2':0.,'removed_P2':0.,
                         'new_only_payload_pairs':0,'removed_payload_pairs':0,'common_payload_pairs':0,
                         'example_new_only_shot':None}
        g=groups[key];fb,fa=bool(failed[0][shot]),bool(failed[1][shot]);weight=float(w[shot])
        g['baseline_P2']+=weight*fb;g['alternating_P2']+=weight*fa
        if fa and not fb:
            g['new_only_P2']+=weight;g['new_only_payload_pairs']+=1
            if g['example_new_only_shot'] is None:g['example_new_only_shot']=int(shot)
        elif fb and not fa:g['removed_P2']+=weight;g['removed_payload_pairs']+=1
        else:g['common_payload_pairs']+=1
        types[typepair]+=np.array([fb,fa])*weight
        for i,ls in enumerate([(la,lb),(al,bl)]):
            ls=sorted(ls,key=lambda l:(l[0]=='B',l[1]))
            position_key=(typepair,ls[0][0]+str(ls[0][1]),ls[1][0]+str(ls[1][1]))
            positions[position_key][i]+=bool(failed[i][shot])*weight
    classes=[]
    for g in groups.values():
        g['delta_P2']=g['alternating_P2']-g['baseline_P2']
        g['delta_P2_survival_corrected']=g['delta_P2']*survival
        classes.append(g)
    classes.sort(key=lambda g:g['delta_P2'],reverse=True)
    def example(shot):
        a,b=int(aa[shot]),int(bb[shot]);payloads=[];replays=[]
        for i,sampler in enumerate(samplers):
            ls=[a,b] if i==0 else [int(mapping[a]),int(mapping[b])]
            faults={};records=[]
            for li,code in zip(ls,[va[shot],vb[shot]]):
                loc=sampler.locations[li];payload=list(sampler.channels[li])[int(code)-1]
                faults.setdefault(loc[:2],{})[loc[2]]=payload
                records.append({'physical_identity':identity(designs[i],loc),'executed_location':loc,'payload':payload})
            replay=detailed_replay(designs[i],faults)
            assert replay['logical_failure']==[bool(traces[i]['x_fail'][shot]),bool(traces[i]['z_fail'][shot])]
            payloads.append(records);replays.append(replay)
        return {'weight_P2':float(w[shot]),'survival_weight':float(w[shot]*survival),
                'baseline_faults':payloads[0],'alternating_faults':payloads[1],
                'baseline_replay':replays[0],'alternating_replay':replays[1]}
    examples=[];seen_types=set()
    for rank,g in enumerate(classes):
        if g['delta_P2']<=0 or g['example_new_only_shot'] is None:continue
        # Top 5 classes plus one representative of each additional broad type.
        if len(examples)<5 or g['location_types'] not in seen_types:
            examples.append({'class_rank':rank+1,'class':g,'fault_propagation':example(g['example_new_only_shot'])})
            seen_types.add(g['location_types'])
    raw=[float(w[f].sum()) for f in failed]
    approved=json.loads((ROOT/'results/stageA_check_order.json').read_text())['rows']
    for name,value in zip(('ZZZXXX','ZXZXZX'),raw):
        assert np.isclose(value,approved['c1_'+name]['counts']['P2'],rtol=3e-12,atol=0)
    delta=raw[1]-raw[0]
    assert np.isclose(partition['new_only']['P2_weight']-partition['removed']['P2_weight'],delta,rtol=3e-12,atol=0)
    out={'stage':'A3','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'preset':f'UniformNoise({p:g})','p':p,'pair_payload_trials_per_design':n,
         'identity_convention':'Same stage, check kind/support, original op index and Pauli payload; executed positions are remapped.',
         'baseline_P2':raw[0],'alternating_P2':raw[1],'delta_P2':delta,
         'survival_factor':survival,'baseline_P2_survival_corrected':raw[0]*survival,
         'alternating_P2_survival_corrected':raw[1]*survival,'delta_P2_survival_corrected':delta*survival,
         'matched_payload_partition':partition,
         'by_location_types':[{'types':key,'baseline_P2':float(v[0]),'alternating_P2':float(v[1]),'delta_P2':float(v[1]-v[0])} for key,v in sorted(types.items())],
         'by_executed_positions':[{'types':key[0],'first_position':key[1],'second_position':key[2],
                                  'baseline_P2':float(v[0]),'alternating_P2':float(v[1]),'delta_P2':float(v[1]-v[0])}
                                 for key,v in sorted(positions.items(),key=lambda item:item[1][1]-item[1][0],reverse=True)],
         'classes_positive_delta':[g for g in classes if g['delta_P2']>1e-18],
         'classes_negative_delta':[g for g in classes if g['delta_P2']< -1e-18],
         'equal_class_count':sum(abs(g['delta_P2'])<=1e-18 for g in classes),
         'top_examples':examples,'scalar_replays_checked':2*len(examples),
         'one_round_MC_comparison':{name:{'p_L':approved['c1_'+name]['monte_carlo']['p_L'],
             'raw_P2_over_MC':approved['c1_'+name]['counts']['P2']/approved['c1_'+name]['monte_carlo']['p_L'],
             'survival_P2_over_MC':approved['c1_'+name]['counts']['P2_survival_corrected']/approved['c1_'+name]['monte_carlo']['p_L']}
             for name in ('ZZZXXX','ZXZXZX')}}
    save(ROOT/'results/stageA3_pair_mechanism.json',out)
    print(json.dumps({k:out[k] for k in ('pair_payload_trials_per_design','baseline_P2','alternating_P2','delta_P2','matched_payload_partition','by_location_types')},indent=2),flush=True)
    return out


if __name__=='__main__':run_mechanism()
