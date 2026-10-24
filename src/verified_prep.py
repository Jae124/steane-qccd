"""B1 only: published Goto ancilla preparation and its Hadamard dual.

Goto, Sci. Rep. 6, 19578 (2016), Fig. 1(c), DOI:10.1038/srep19578.
Paper local lines 0..6 use the same code as steane_ft.SUPPORTS (same row span).
Default physical layout reserves data 0..6, block 7..13 and verifier 14.
Only the 8-qubit preparation gadget is run here; no EC coupling is implemented.
"""
from dataclasses import dataclass
import itertools
import math
import numpy as np
import stim

from flag_design import layer_timing
from ion_noise import IonNoise,IdleIonNoise
from ion_timing import idle_probability
import steane_ft as ft
import steane_noise as scalar

BLOCK=tuple(range(7,14));VERIFIER=14;N_QUBITS=15
PAPER_PLUS_INPUTS=(1,2,3)
PAPER_ENCODER=((1,0),(3,5),(2,6),(1,4),(2,0),(3,6),(1,5),(6,4))
VERIFY_SUPPORT=(0,5,6)
TYPES=('reset','measurement','2q','idle')
STABILIZERS=tuple(sum(1<<q for q in support) for support in ft.SUPPORTS)
CODE_GROUP=tuple((STABILIZERS[0] if i&1 else 0)^(STABILIZERS[1] if i&2 else 0)^
                 (STABILIZERS[2] if i&4 else 0) for i in range(8))
PREP_PHASE_GROUP=tuple(sorted(set(CODE_GROUP)|{s^127 for s in CODE_GROUP}))
def popcount(value):return bin(int(value)).count('1')

REDUCED_WEIGHT=np.array([min(popcount(v^s) for s in CODE_GROUP) for v in range(128)],np.uint8)
STATE_PHASE_WEIGHT=np.array([min(popcount(v^s) for s in PREP_PHASE_GROUP) for v in range(128)],np.uint8)
RAW_WEIGHT=np.array([popcount(v) for v in range(128)],np.uint8)


def preparation_ops(basis,*,block=BLOCK,verifier=VERIFIER):
    """Fixed Figure 1(c) transcription; + is its full H-conjugate, not a new fit.

    Disjoint CNOTs drawn in the same paper column are linearized top-to-bottom.
    All gates remain serial in the timing model. All 8 initial resets precede
    the encoder; this convention is fixed before single-fault validation.
    """
    if basis not in ('0','+'):raise ValueError('basis must be 0 or +')
    if len(block)!=7 or len(set((*block,verifier)))!=8:raise ValueError('Use 7 distinct block qubits and a distinct verifier')
    ops=[('RX' if q in PAPER_PLUS_INPUTS else 'R',block[q]) for q in range(7)]
    ops.append(('R',verifier))
    ops += [('CX',block[c],block[t]) for c,t in PAPER_ENCODER]
    ops += [('CX',block[q],verifier) for q in VERIFY_SUPPORT]
    ops.append(('M',verifier,'verify'))
    if basis=='+':
        dual=[]
        for op in ops:
            if op[0]=='CX':dual.append(('CX',op[2],op[1]))
            elif op[0] in ('R','RX'):dual.append(('RX' if op[0]=='R' else 'R',op[1]))
            else:dual.append(('MX',op[1],op[2]))
        ops=dual
    return ops


def op_part(index):
    return 'initial_resets' if index<8 else 'encoder' if index<16 else 'verification' if index<19 else 'verification_readout'


def insert_preparation_idles(ops,params,policy,*,output_qubits=BLOCK,live_on_entry=()):
    """Existing SEQ/PAR_MR rule with preparation live windows and output waits.

    Standalone B1 has no data live on entry. A later embedding must explicitly
    include live spectator data and charge waiting during every retry. Waits
    split at prep-part boundaries without adding physical time.
    """
    starts,ends,total,readout=layer_timing(ops,params,policy)
    boundaries=sorted(set([0.,starts[8],starts[16],starts[19],total]))
    last={q:0. for q in live_on_entry};out=[];metadata=[];mapping={}
    def idle(q,start,end,before):
        cuts=[start]+[v for v in boundaries if start+1e-15<v<end-1e-15]+[end]
        for lo,hi in zip(cuts,cuts[1:]):
            if hi-lo<=1e-15:continue
            mid=(lo+hi)/2
            part=('initial_resets' if mid<starts[8] else 'encoder' if mid<starts[16]
                  else 'verification' if mid<starts[19] else 'verification_readout')
            out.append(('IDLE',q,hi-lo));metadata.append({'part':part,'original_op_index':None,'before_original_op':before})
    for idx,op in enumerate(ops):
        touched=op[1:3] if op[0]=='CX' else (op[1],)
        for q in touched:
            if q in last:idle(q,last[q],starts[idx],idx)
        mapping[idx]=len(out);out.append(op);metadata.append({'part':op_part(idx),'original_op_index':idx})
        for q in touched:
            last[q]=ends[idx]
            if op[0] in ('M','MX'):last.pop(q,None)
    for q in sorted(set(output_qubits)|set(live_on_entry)):idle(q,last[q],total,len(ops))
    return out,metadata,mapping,total,readout


def _pauli(sim,qubits,payload):
    for q,p in zip(qubits,payload):
        if p=='X':sim.x(q)
        elif p=='Y':sim.y(q)
        elif p=='Z':sim.z(q)


def run_tableau(ops,faults=None,*,sim=None):
    """New dynamic-15-qubit primitive; the frozen 9-qubit executor is untouched."""
    sim=stim.TableauSimulator() if sim is None else sim
    faults=faults or {};bits={}
    for idx,op in enumerate(ops):
        kind,q=op[:2];payload=faults.get(idx)
        if kind=='R':sim.reset(q)
        elif kind=='RX':sim.reset_x(q)
        elif kind=='CX':sim.cnot(q,op[2])
        elif kind=='IDLE':pass
        elif kind in ('M','MX'):
            if kind=='MX':sim.h(q)
            bit=sim.measure(q)
            if kind=='MX':sim.h(q)
            bits[op[2]]=int(bit)^(payload=='FLIP')
            continue
        else:raise ValueError(f'Unknown prep operation {op}')
        if payload is not None:_pauli(sim,op[1:3] if kind=='CX' else [q],payload)
    return sim,bits


def _observable(support,pauli):
    result=stim.PauliString(N_QUBITS)
    for q in support:result[q]=pauli
    return result


def tableau_signature(sim,basis):
    return tuple(sim.peek_observable_expectation(_observable([BLOCK[q] for q in support],p))
                 for p in ('Z','X') for support in ft.SUPPORTS)+(sim.peek_observable_expectation(_observable(BLOCK,'Z' if basis=='0' else 'X')),)


def frame_signature(fx,fz,basis):
    x,z=(fx>>7)&127,(fz>>7)&127
    return tuple(1-2*(popcount(v&s)&1) for v in (x,z) for s in STABILIZERS)+(1-2*(popcount(x if basis=='0' else z)&1),)


def weights(fx,fz,basis):
    x,z=(fx>>7)&127,(fz>>7)&127
    danger,other=(x,z) if basis=='0' else (z,x)
    return {'dangerous_code_reduced':int(REDUCED_WEIGHT[danger]),
            'other_raw':int(RAW_WEIGHT[other]),'other_code_reduced':int(REDUCED_WEIGHT[other]),
            'other_prepared_state_reduced':int(STATE_PHASE_WEIGHT[other]),
            'dangerous_mask':danger,'other_mask':other}


def dangerous_weight_from_signature(signature,basis):
    bits=tuple(v==-1 for v in ((signature[:3] if basis=='0' else signature[3:6])+(signature[6],)))
    return min(popcount(mask) for mask in range(128)
               if tuple(bool(popcount(mask&s)&1) for s in STABILIZERS)+(bool(popcount(mask)&1),)==bits)


class VerifiedPreparation:
    def __init__(self,basis,noise,*,timing_policy='SEQ'):
        self.basis=basis;self.noise=noise;self.params=getattr(noise,'params',None)
        self.timed=isinstance(noise,IdleIonNoise)
        if getattr(noise,'transport_events',()):raise ValueError('B1 transport stays empty')
        if timing_policy not in ('SEQ','PAR_MR'):raise ValueError('Use SEQ or PAR_MR')
        if not self.timed and timing_policy!='SEQ':raise ValueError('Untimed presets have no PAR_MR timing')
        self.policy=timing_policy if self.timed else 'timing-free'
        self.base_ops=preparation_ops(basis)
        if self.timed:
            self.ops,self.metadata,self.mapping,self.duration,self.readout=insert_preparation_idles(self.base_ops,self.params,timing_policy)
        else:
            self.ops=self.base_ops
            self.metadata=[{'part':op_part(i),'original_op_index':i} for i in range(len(self.ops))]
            self.mapping={i:i for i in range(len(self.ops))};self.duration=self.readout=None
        self.locations=[('P',0,i,op) for i,op in enumerate(self.ops)]
        self.channels=[scalar.validated_channel(self,loc) for loc in self.locations]
        self.influences=[]
        for i,ch in enumerate(self.channels):
            values=[]
            for payload in ch:
                fx,fz,bits=scalar.run_ops_frame(0,0,self.ops,{i:payload})
                danger=((fx if basis=='0' else fz)>>7)&127
                values.append({'payload':payload,'fx':fx,'fz':fz,'verify':bits['verify'],
                               'code':danger|(bits['verify']<<7)})
            self.influences.append(values)

    def channel(self,loc):
        if loc[3][0]=='IDLE':
            p=idle_probability(loc[3][2],self.params.T2)
            return {'Z':p} if p else {}
        return IonNoise.channel(self.noise,loc) if self.timed else self.noise.channel(loc)

    def run_frame(self,faults=None,*,fx=0,fz=0):
        fx,fz,bits=scalar.run_ops_frame(fx,fz,self.ops,faults or {})
        return fx,fz,bits['verify']

    def noiseless_check(self):
        sim,bits=run_tableau(self.ops)
        fx,fz,verify=self.run_frame()
        assert bits['verify']==verify==0
        assert tableau_signature(sim,self.basis)==frame_signature(fx,fz,self.basis)==(1,)*7
        return {'verification_outcome':0,'signature':[1]*7}

    def check_faults(self,*,universal=False):
        """No circuit/decoder changes are allowed after a genuine FT failure."""
        failures=[];mismatches=[];runs=accepted=rejected=0;parts={};maxima={'other_raw':0,'other_code_reduced':0,'other_prepared_state_reduced':0,'dangerous_code_reduced':0}
        for idx,op in enumerate(self.ops):
            choices=([payload for _,payload in ft.fault_list([op])] if universal else
                     [v for v,p in self.channels[idx].items() if p>0])
            for payload in choices:
                fx,fz,bit=self.run_frame({idx:payload})
                sim,bits=run_tableau(self.ops,{idx:payload})
                signature=frame_signature(fx,fz,self.basis);tableau=tableau_signature(sim,self.basis)
                w=weights(fx,fz,self.basis);tw=dangerous_weight_from_signature(tableau,self.basis)
                record={'gadget':'verified_'+self.basis,'part':self.metadata[idx]['part'],'op_index':idx,
                        'original_op_index':self.metadata[idx]['original_op_index'],'op':op,'payload':payload,
                        'verification_outcome':bit,'accepted':bit==0,'syndromes_and_logical_signature':signature,
                        'correction_applied':[],'residual_frame':{'X':fx,'Z':fz},'weights':w}
                if bit!=bits['verify'] or signature!=tableau or w['dangerous_code_reduced']!=tw:
                    mismatches.append(dict(record,tableau_verification=bits['verify'],tableau_signature=tableau,tableau_dangerous_weight=tw))
                if bit==0:
                    accepted+=1
                    for key in maxima:maxima[key]=max(maxima[key],w[key])
                    if w['dangerous_code_reduced']>1:failures.append(record)
                else:rejected+=1
                runs+=1;part=parts.setdefault(self.metadata[idx]['part'],{'runs':0,'accepted':0,'rejected':0})
                part['runs']+=1;part['accepted']+=bit==0;part['rejected']+=bit!=0
        return {'catalog':'universal physical reset/CX/measurement/idle faults' if universal else 'positive preset-channel payloads',
                'prep_fault_runs':runs,'accepted_single_fault_runs':accepted,'rejected_single_fault_runs':rejected,
                'accepted_dangerous_failures':len(failures),'frame_tableau_mismatches':len(mismatches),
                'accepted_maxima':maxima,'by_part':parts,'failure_details':failures,'mismatch_details':mismatches}

    def compiled_codes(self):
        return [np.array([0]+[v['code'] for v in table],np.uint8) for table in self.influences]

    def exact_distribution(self):
        """Exact finite-rate XOR convolution over danger-mask + verification bit."""
        prob=np.zeros(256);prob[0]=1.;states=np.arange(256)
        for ch,table in zip(self.channels,self.influences):
            previous=prob;prob=(1-sum(ch.values()))*previous
            for p,v in zip(ch.values(),table):prob=prob+p*previous[states^v['code']]
        assert abs(prob.sum()-1)<1e-12
        accepted=np.arange(256)<128;bad=accepted&(REDUCED_WEIGHT[np.arange(256)&127]>1)
        acceptance=float(prob[accepted].sum());hazard=float(prob[bad].sum())
        return {'acceptance':acceptance,'rejection':1-acceptance,'accepted_dangerous_probability':hazard,
                'dangerous_given_acceptance':hazard/acceptance,'method':'Exact independent-location Pauli XOR convolution; no truncation.'}

    def count_pairs(self):
        raw=corrected=first=0.;by_type={};by_part={};trials=malignant=0
        rates=[sum(ch.values()) for ch in self.channels]
        for ch,table in zip(self.channels,self.influences):
            for p,v in zip(ch.values(),table):
                if not v['verify'] and REDUCED_WEIGHT[v['code']&127]>1:first+=p
        for a,b in itertools.combinations(range(len(self.ops)),2):
            survival=math.prod(1-rate for i,rate in enumerate(rates) if i not in (a,b))
            kind=lambda i:'measurement' if self.ops[i][0] in ('M','MX') else 'reset' if self.ops[i][0] in ('R','RX') else '2q' if self.ops[i][0]=='CX' else 'idle'
            typepair='*'.join(sorted([kind(a),kind(b)]));partpair='*'.join(sorted([self.metadata[a]['part'],self.metadata[b]['part']]))
            for pa,va in zip(self.channels[a].values(),self.influences[a]):
                for pb,vb in zip(self.channels[b].values(),self.influences[b]):
                    if pa<=0 or pb<=0:continue
                    trials+=1;code=va['code']^vb['code']
                    if code<128 and REDUCED_WEIGHT[code&127]>1:
                        value=pa*pb;raw+=value;corrected+=value*survival;malignant+=1
                        by_type[typepair]=by_type.get(typepair,0.)+value
                        by_part[partpair]=by_part.get(partpair,0.)+value
        mean=sum(rates)/len(rates)
        return {'definition':'Per-attempt probability mass of accepted preparations with dangerous component reduced weight >1; not full-EC p_L.',
                'P1':first,'P2':raw,'P2_survival_corrected':corrected,
                'P2_survival_mean_approx':raw*(1-mean)**(len(rates)-2),
                'positive_payload_pair_trials':trials,'malignant_payload_pairs':malignant,
                'P2_by_location_type_pair':by_type,'P2_by_prep_part_pair':by_part,
                'n_locations':len(rates),'mean_location_fault_probability':mean}

    def sample_attempts(self,shots,rng,*,return_arrays=False):
        """Exact finite-rate sparse location draws; no quantum decisions use raw codewords."""
        codes=np.zeros(shots,np.uint8);signals=np.zeros((4,shots),np.uint8)
        fx=np.zeros(shots,np.uint16) if return_arrays else None
        fz=np.zeros(shots,np.uint16) if return_arrays else None
        compiled=self.compiled_codes()
        for idx,(ch,lookup) in enumerate(zip(self.channels,compiled)):
            probabilities=np.array(list(ch.values()),float);p=float(probabilities.sum())
            hits=rng.binomial(shots,p)
            if not hits:continue
            indices=rng.choice(shots,hits,replace=False)
            payload=rng.choice(len(probabilities),hits,p=probabilities/p)+1
            codes[indices]^=lookup[payload]
            if return_arrays:
                fx[indices]^=np.array([0]+[v['fx'] for v in self.influences[idx]],np.uint16)[payload]
                fz[indices]^=np.array([0]+[v['fz'] for v in self.influences[idx]],np.uint16)[payload]
            op=self.ops[idx][0]
            typ='measurement' if op in ('M','MX') else 'reset' if op in ('R','RX') else '2q' if op=='CX' else 'idle'
            signals[TYPES.index(typ),indices]^=lookup[payload]>>7
        reject=(codes>>7).astype(bool);bad=(~reject)&(REDUCED_WEIGHT[codes&127]>1)
        count=signals.sum(axis=0)
        causes={t:int(np.count_nonzero(reject&(count==1)&signals[i].astype(bool))) for i,t in enumerate(TYPES)}
        causes['mixed']=int(np.count_nonzero(reject&(count>1)))
        out={'shots':int(shots),'accepted':int((~reject).sum()),'rejected':int(reject.sum()),
             'accepted_dangerous':int(bad.sum()),'verification_rejection_cause_counts':causes,
             'verification_contributing_type_counts':{t:int(np.count_nonzero(reject&signals[i].astype(bool))) for i,t in enumerate(TYPES)},
             'verification_fractional_type_counts':{t:float(np.sum((reject&signals[i].astype(bool))/np.maximum(count,1))) for i,t in enumerate(TYPES)}}
        if return_arrays:return out,{'codes':codes,'fx':fx,'fz':fz,'rejected':reject,'accepted_dangerous':bad}
        return out

    def sample_retries(self,shots,rng,*,max_attempts):
        """Explicit total-attempt cap; no interpretation of max_retries is hidden."""
        if type(max_attempts) is not int or max_attempts<1:raise ValueError('Use a positive total-attempt cap')
        active=shots;hist=[];hazards=executions=0
        for _ in range(max_attempts):
            executions+=active
            if active:
                counts=self.sample_attempts(active,rng)
                hist.append(counts['accepted']);hazards+=counts['accepted_dangerous'];active=counts['rejected']
            else:hist.append(0)
        return {'shots':shots,'max_attempts':max_attempts,'accepted_on_attempt':hist,'all_attempts_rejected':active,
                'accepted_dangerous':hazards,'standalone_failures':active+hazards,'attempt_executions':executions,
                'mean_attempts':executions/shots,
                'scope':'Standalone preparations only; no data block idle or EC coupling has been simulated.'}
