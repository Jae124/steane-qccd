"""Steane-style adapter for the approved A3 memory driver; B2 smoke validation.

Whole Clifford segment influences are compiled from the independent scalar
frame executor. Incoming data survive retries/rounds. Aborted shots are
absorbing and execute no more operations. No fit/comparison is performed here.
"""
import numpy as np
import steane_ft as ft
import steane_noise as scalar
from memory_mc import RoundStep,flag_endpoint_score
from steane_extraction import syndrome

SYNDROME=np.array([syndrome([(v>>q)&1 for q in range(7)]) for v in range(128)],np.uint8)
FIX=np.array([sum(1<<q for q in ft.standard_decode(s)) for s in range(8)],np.uint16)
TYPES=('reset','measurement','2q','idle')


class SteaneMemoryRound:
    branch_label='verification_rejection'
    def __init__(self,design):
        self.design=design;self.compiled={}
        for key,seg in design.segments.items():
            table=[]
            for i,op in enumerate(seg['ops']):
                channel=design.channel((key,i,op));variants=list(channel);xs=[0];zs=[0];signals=[0]
                for v in variants:
                    x,z,bits=scalar.run_ops_frame(0,0,seg['ops'],{i:v})
                    signal=bits['verify'] if key[0]=='P' else syndrome(bits)
                    xs.append(x);zs.append(z);signals.append(signal)
                typ='reset' if op[0] in ('R','RX') else '2q' if op[0]=='CX' else 'measurement' if op[0] in ('M','MX') else 'idle'
                table.append({'variants':variants,'probabilities':np.array(list(channel.values())),
                    'xs':np.array(xs,np.uint16),'zs':np.array(zs,np.uint16),'signals':np.array(signals,np.uint8),'type':TYPES.index(typ)})
            self.compiled[key]=table

    def initial_state(self,shots):
        return np.zeros(shots,np.uint16),np.zeros(shots,np.uint16),np.ones(shots,bool)

    def explicit_faults(self,shot_faults):
        """Scalar payload dictionaries -> sparse codes, for exact replay tests."""
        out={}
        for shot,faults in enumerate(shot_faults):
            for key,entries in faults.items():
                for i,payload in entries.items():
                    hits,codes=out.setdefault(key,{}).setdefault(i,([],[]))
                    hits.append(shot);codes.append(self.compiled[key][i]['variants'].index(payload)+1)
        return {key:{i:(np.array(h,int),np.array(c,np.uint8)) for i,(h,c) in entries.items()} for key,entries in out.items()}

    def step(self,state,rng,*,faults=None):
        x,z,alive=(a.copy() for a in state);shots=len(x)
        branch=np.zeros(shots,bool);aborted=np.zeros(shots,bool);reject_causes={t:0 for t in (*TYPES,'mixed')}
        executions={str(k):0 for k in self.design.segments};rejections=attempts=0
        def segment(key,active):
            n=len(active);signal=np.zeros(n,np.uint8);causes=np.zeros((4,n),np.uint8)
            if key[0]=='P':x[active]&=127;z[active]&=127
            elif self.design.order[key[1]]=='X':
                signal[:]=SYNDROME[((x[active]>>7)^x[active])&127]
                z[active]^=(z[active]>>7)&127
            else:
                signal[:]=SYNDROME[((z[active]>>7)^z[active])&127]
                x[active]^=(x[active]>>7)&127
            executions[str(key)]+=n
            for i,t in enumerate(self.compiled[key]):
                if faults is None:
                    p=float(t['probabilities'].sum());nh=int(rng.binomial(n,min(p,1.)))
                    if not nh:continue
                    local=rng.choice(n,nh,replace=False)
                    code=rng.choice(len(t['variants']),nh,p=t['probabilities']/p)+1
                else:
                    values=faults.get(key,{}).get(i)
                    if values is None:continue
                    hits,codes=values;pos=np.searchsorted(active,hits);ok=pos<n
                    ok[ok]&=active[pos[ok]]==hits[ok]
                    local=pos[ok];code=codes[ok]
                selected=active[local];x[selected]^=t['xs'][code];z[selected]^=t['zs'][code]
                signal[local]^=t['signals'][code];causes[t['type'],local]^=t['signals'][code]
            if key[0]=='C':x[active]&=127;z[active]&=127
            return signal,causes
        for hi,kind in enumerate(self.design.order):
            waiting=np.flatnonzero(alive);accepted=[]
            for a in range(self.design.max_attempts):
                if not len(waiting):break
                signal,causes=segment(('P',hi,a),waiting);attempts+=len(waiting)
                rejected=signal.astype(bool);branch[waiting[rejected]]=True;rejections+=int(rejected.sum())
                multiplicity=(causes[:,rejected]!=0).sum(axis=0)
                for ti,t in enumerate(TYPES):reject_causes[t]+=int(np.count_nonzero((causes[ti,rejected]!=0)&(multiplicity==1)))
                reject_causes['mixed']+=int(np.count_nonzero(multiplicity>1))
                accepted.extend(waiting[~rejected]);waiting=waiting[rejected]
            alive[waiting]=False;aborted[waiting]=True
            active=np.sort(np.array(accepted,int))
            if len(active):
                signal,_=segment(('C',hi,0),active)
                if kind=='X':x[active]^=FIX[signal]
                else:z[active]^=FIX[signal]
        return RoundStep((x&127,z&127,alive),branch,aborted,
            {'attempts':attempts,'rejections':rejections,'new_aborts':int(aborted.sum()),
             'segment_execution_counts':executions,'rejection_cause_counts':reject_causes})


def endpoint_score(state):
    return flag_endpoint_score(state[:2])
