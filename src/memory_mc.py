"""Reusable noisy-round memory harness: carry state; score only at the endpoint.

An adapter provides initial_state(shots) and step(state, rng, faults=None) ->
RoundStep. The driver does not know the circuit, decoder, ancilla count or retry
policy. A future adapter can report an aborted preparation as an absorbing shot
failure. FlagMemoryRound is the only physical adapter implemented in A3.
"""
from dataclasses import dataclass
from typing import Protocol
import math
import numpy as np

from design_mc import BatchSampler, LOGICAL

CATEGORIES=('reset','measurement','2q','idle','carry')


@dataclass
class RoundStep:
    state: object
    branch_taken: np.ndarray
    aborted: np.ndarray
    diagnostics: dict


class RoundAdapter(Protocol):
    def initial_state(self, shots): ...
    def step(self, state, rng, *, faults=None): ...


class FlagMemoryRound:
    """Existing flag circuit/channels/decoder, accepting an arbitrary input frame."""
    branch_label='bare'
    def __init__(self,design):
        self.design=design
        self.sampler=BatchSampler(design)

    def initial_state(self,shots):
        return np.zeros(shots,np.uint16),np.zeros(shots,np.uint16)

    def step(self,state,rng,*,faults=None):
        incoming_x,incoming_z=state
        shots=len(incoming_x)
        if faults is None:faults=self.sampler.sample_faults(shots,rng)
        pieces=[np.flatnonzero(incoming_x | incoming_z)] + [v[0] for v in faults.values()]
        indices=np.unique(np.concatenate(pieces));n=len(indices)
        codes=np.zeros((len(self.sampler.locations),n),np.uint8)
        for li,(hits,payloads) in faults.items():
            codes[li,np.searchsorted(indices,hits)]=payloads
        fx=np.zeros((5,n),np.uint16);fz=np.zeros_like(fx)
        fx[4]=incoming_x[indices];fz[4]=incoming_z[indices]
        stop=np.full(n,-1,np.int8);flag=np.zeros(n,np.uint8);trigger_synd=np.zeros(n,np.uint8)
        contributors=np.zeros((5,n),bool)

        def run(key,active):
            x,z=fx[:,active].copy(),fz[:,active].copy();bits={}
            for li in self.sampler.by_circuit[key]:
                op=self.sampler.locations[li][3];kind,q=op[:2]
                category,xs,zs,flips=self.sampler.compiled[li];code=codes[li,active]
                if kind in ('R','RX'):
                    x &= np.uint16(65535^(1<<q));z &= np.uint16(65535^(1<<q))
                elif kind=='CX':
                    t=op[2];x ^= ((x>>q)&1)<<t;z ^= ((z>>t)&1)<<q
                elif kind in ('M','MX'):
                    bit=((x if kind=='M' else z)>>q)&1
                    bit[category] ^= flips[code];bits[op[2]]=bit
                    continue
                x[category] ^= xs[code];z[category] ^= zs[code]
            fx[:,active],fz[:,active]=x,z
            return bits

        for i in range(6):
            active=np.flatnonzero(stop<0)
            if not len(active):break
            bits=run(('F',i),active)
            synd=np.bitwise_xor.reduce(bits['synd'],axis=0)
            flagbit=np.bitwise_xor.reduce(bits['flag'],axis=0)
            hit=(synd|flagbit).astype(bool);chosen=active[hit]
            stop[chosen]=i;flag[chosen]=flagbit[hit];trigger_synd[chosen]=synd[hit]
            contributors[:,chosen]=((bits['synd'][:,hit]&synd[hit])|(bits['flag'][:,hit]&flagbit[hit])).astype(bool)
        active=np.flatnonzero(stop>=0);sx=np.zeros(n,np.uint8);sz=np.zeros(n,np.uint8)
        for i,(kind,k) in enumerate(self.design.round_order):
            if not len(active):break
            bits=run(('B',i),active)
            syndrome=np.bitwise_xor.reduce(bits['synd'],axis=0).astype(np.uint8)
            (sx if kind=='Z' else sz)[active] |= syndrome<<(2-k)
        x=np.bitwise_xor.reduce(fx,axis=0);z=np.bitwise_xor.reduce(fz,axis=0)
        xfix=self.sampler.standard[sx];zfix=self.sampler.standard[sz]
        for i,(kind,k) in enumerate(self.design.round_order):
            chosen=(stop==i)&(flag!=0)
            if kind=='Z':zfix[chosen]=self.sampler.flagged[k][sz[chosen]]
            else:xfix[chosen]=self.sampler.flagged[k][sx[chosen]]
        x ^= xfix;z ^= zfix
        out_x=incoming_x.copy();out_z=incoming_z.copy()
        out_x[indices]=x;out_z[indices]=z
        branch=np.zeros(shots,bool);branch[indices]=stop>=0
        multiplicity=contributors.sum(axis=0)
        if np.any((stop>=0)&(multiplicity==0)):raise RuntimeError('Trigger has no signal contributor')
        exclusive={c:int(np.count_nonzero(contributors[i]&(multiplicity==1))) for i,c in enumerate(CATEGORIES)}
        exclusive['mixed']=int(np.count_nonzero(multiplicity>1))
        diag={'triggers':int(branch.sum()),'trigger_circuit_counts':[int(np.count_nonzero(stop==i)) for i in range(6)],
              'trigger_cause_counts':exclusive,
              'trigger_contributing_type_counts':{c:int(contributors[i].sum()) for i,c in enumerate(CATEGORIES)},
              'trigger_fractional_type_counts':{c:float(np.sum(contributors[i]/np.maximum(multiplicity,1))) for i,c in enumerate(CATEGORIES)},
              'trace':{'shot_indices':indices,'stop':stop,'flag':flag,'trigger_synd':trigger_synd,
                       's_x':sx,'s_z':sz,'x_fix':xfix,'z_fix':zfix,'contributors':contributors}}
        return RoundStep((out_x,out_z),branch,np.zeros(shots,bool),diag)


def flag_endpoint_score(state):
    """Ideal decode for observation only: no mutation or intermediate correction."""
    x,z=state
    return LOGICAL[x&127].astype(bool),LOGICAL[z&127].astype(bool)


def _add_diagnostics(total,addition):
    if not total:
        return {k:v.copy() if isinstance(v,(dict,list)) else v for k,v in addition.items() if k!='trace'}
    for key,value in addition.items():
        if key=='trace':continue
        if isinstance(value,dict):
            for label,count in value.items():total[key][label]+=count
        elif isinstance(value,list):total[key]=[a+b for a,b in zip(total[key],value)]
        else:total[key]+=value
    return total


def memory_trial(adapter:RoundAdapter,shots,rounds,rng,*,score=flag_endpoint_score,fault_rounds=None,trace=False):
    """No reset, ideal correction, scoring or extra time between round steps.

    Aborts are absorbing endpoint failures (for future verified-prep adapters).
    Independent R-point streams are used by the experiment's binomial fits.
    """
    if shots<=0 or rounds<=0:raise ValueError('Positive shots and rounds required')
    state=adapter.initial_state(shots);any_branch=np.zeros(shots,bool);aborted=np.zeros(shots,bool)
    diagnostics={};round_diagnostics=[];traces=[]
    for r in range(rounds):
        step=adapter.step(state,rng,faults=None if fault_rounds is None else fault_rounds[r])
        state=step.state;any_branch |= step.branch_taken;aborted |= step.aborted
        diagnostics=_add_diagnostics(diagnostics,step.diagnostics)
        round_diagnostics.append({k:v for k,v in step.diagnostics.items() if k!='trace'})
        if trace:traces.append({'state':state,'branch_taken':step.branch_taken,'diagnostics':step.diagnostics})
    xfail,zfail=score(state)
    fail_any=xfail|zfail|aborted
    union=int(fail_any.sum());xcount=int(np.count_nonzero(xfail|aborted));zcount=int(np.count_nonzero(zfail|aborted))
    counts={'shots':int(shots),'rounds':int(rounds),'fails':union,
            'fails_X':xcount,'fails_Z':zcount,'aborted':int(aborted.sum()),
            'logical_failure_counts':{'X_only':int(np.count_nonzero(xfail&~zfail&~aborted)),
                                      'Z_only':int(np.count_nonzero(zfail&~xfail&~aborted)),
                                      'both_XZ':int(np.count_nonzero(xfail&zfail&~aborted))},
            'any_branch_shots':int(any_branch.sum()),
            'failures_through_any_branch':int(np.count_nonzero(fail_any&any_branch)),
            'failures_through_last_branch':int(np.count_nonzero(fail_any&step.branch_taken)),
            'round_totals':diagnostics,'per_round':round_diagnostics}
    if trace:return counts,{'final_state':state,'rounds':traces,'fail_X':xfail,'fail_Z':zfail}
    return counts


def merge_memory_counts(total,addition):
    if not total:return addition
    if total['rounds']!=addition['rounds']:raise ValueError('Cannot merge different memory lengths')
    for key,value in addition.items():
        if key=='rounds':continue
        if key=='round_totals':total[key]=_add_diagnostics(total[key],value)
        elif key=='per_round':
            total[key]=[_add_diagnostics(a,b) for a,b in zip(total[key],value)]
        elif isinstance(value,dict):
            for label,count in value.items():total[key][label]+=count
        else:total[key]+=value
    return total


def memory_metrics(counts):
    n=counts['shots'];denom=n*counts['rounds'];out=dict(counts)
    for suffix,key in [('', 'fails'),('_X','fails_X'),('_Z','fails_Z')]:
        p=counts[key]/n
        out['p_L'+suffix]=p;out['stderr'+suffix]=math.sqrt(max(p,1/n)*(1-p)/n)
    out['basis_labels']={'p_L_X':'|0_L>: logical X parity failure','p_L_Z':'|+_L>: logical Z parity failure',
                         'p_L':'union of either basis failure (secondary if basis fits are primary)'}
    d=counts['round_totals']
    out['trigger_rate_per_round']=d['triggers']/denom
    out['at_least_one_bare_probability']=counts['any_branch_shots']/n
    out['trigger_cause_rates_per_round']={k:v/denom for k,v in d['trigger_cause_counts'].items()}
    out['trigger_contributing_type_rates_per_round']={k:v/denom for k,v in d['trigger_contributing_type_counts'].items()}
    out['trigger_fractional_type_rates_per_round']={k:v/denom for k,v in d['trigger_fractional_type_counts'].items()}
    out['trigger_position_rates_per_round']=[v/denom for v in d['trigger_circuit_counts']]
    out['failure_share_through_any_bare']=counts['failures_through_any_branch']/counts['fails'] if counts['fails'] else None
    out['failure_share_through_last_bare']=counts['failures_through_last_branch']/counts['fails'] if counts['fails'] else None
    return out


def memory_single_fault_checks(adapter,rounds=16,placements=(0,7,15)):
    """Every supported payload alone, at early/middle/last round positions."""
    injections=[{} for _ in range(rounds)];shot=0;per_placement={};stage_runs={'F':0,'B':0}
    for r in placements:
        start=shot
        for li,ch in enumerate(adapter.sampler.channels):
            for code,prob in enumerate(ch.values(),1):
                if prob<=0:continue
                hits,codes=injections[r].setdefault(li,([],[]))
                hits.append(shot);codes.append(code);shot+=1
                stage_runs[adapter.sampler.locations[li][0]]+=1
        per_placement[str(r+1)]=shot-start
    injections=[{li:(np.array(a),np.array(c)) for li,(a,c) in d.items()} for d in injections]
    counts=memory_trial(adapter,shot,rounds,np.random.default_rng(0),fault_rounds=injections)
    return {'rounds':rounds,'placement_rounds_one_based':[r+1 for r in placements],
            'payload_runs':shot,'both_basis_runs':2*shot,'payload_runs_by_placement':per_placement,
            'flagged_payload_runs':stage_runs['F'],'bare_only_payload_runs':stage_runs['B'],
            'failures_union':counts['fails'],'failures_X':counts['fails_X'],'failures_Z':counts['fails_Z']}
