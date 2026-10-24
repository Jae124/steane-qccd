"""Stage B2: fixed verified ancillas, transversal extraction, standard decoding.

Only H*m and verification parity drive decisions. Ideal codeword readouts are
random; raw measurement bits are deliberately absent from comparison logs.
Rejected preparations reset ancillas, never data. All attempts accrue data idle.
"""
from collections import Counter
import math
import numpy as np
import stim

import steane_ft as ft
import steane_noise as scalar
from flag_design import layer_timing
from ion_noise import IonNoise,IdleIonNoise
from ion_timing import idle_probability
from verified_prep import (BLOCK,VERIFIER,preparation_ops,insert_preparation_idles,
                           run_tableau,popcount)

DATA=tuple(range(7));MASK=127
HALF_ORDERS={'X_FIRST':('X','Z'),'Z_FIRST':('Z','X')}
DEFAULT_MAX_ATTEMPTS=4


def syndrome(bits):
    return ft.syndrome_value([sum(bits[q] for q in support)%2 for support in ft.SUPPORTS])


def frame_signature(x,z,basis):
    return tuple(1-2*(popcount(v&sum(1<<q for q in s))%2)
                 for v in (x,z) for s in ft.SUPPORTS)+(1-2*(popcount((x if basis=='0' else z)&MASK)%2),)


def tableau_signature(sim,basis):
    def observable(s,p):
        o=stim.PauliString(15)
        for q in s:o[q]=p
        return o
    return tuple(sim.peek_observable_expectation(observable(s,p))
                 for p in ('Z','X') for s in ft.SUPPORTS)+(sim.peek_observable_expectation(observable(DATA,'Z' if basis=='0' else 'X')),)


def extraction_ops(kind):
    if kind=='X':return [('CX',q,BLOCK[q]) for q in DATA]+[('M',BLOCK[q],q) for q in DATA]
    if kind=='Z':return [('CX',BLOCK[q],q) for q in DATA]+[('MX',BLOCK[q],q) for q in DATA]
    raise ValueError('Half kind X or Z required')


def extraction_idles(ops,params,policy):
    starts,ends,total,readout=layer_timing(ops,params,policy)
    last={q:0. for q in (*DATA,*BLOCK)};out=[];meta=[]
    boundary=starts[7]
    def wait(q,lo,hi):
        cuts=[lo]+([boundary] if lo+1e-15<boundary<hi-1e-15 else [])+[hi]
        for a,b in zip(cuts,cuts[1:]):
            if b-a>1e-15:
                out.append(('IDLE',q,b-a));meta.append({'part':'coupling' if (a+b)/2<boundary else 'ancilla_readout','original_op_index':None})
    for i,op in enumerate(ops):
        touched=op[1:3] if op[0]=='CX' else (op[1],)
        for q in touched:
            if q in last:wait(q,last[q],starts[i])
        out.append(op);meta.append({'part':'coupling' if i<7 else 'ancilla_readout','original_op_index':i})
        for q in touched:
            last[q]=ends[i]
            if op[0] in ('M','MX'):last.pop(q)
    for q in DATA:wait(q,last[q],total)
    return out,meta,total,readout


class SteaneExtraction:
    """Sequential reusable block; synthetic SEQ/PAR_MR layers, empty transport."""
    def __init__(self,noise,*,half_order='X_FIRST',timing_policy='SEQ',max_attempts=DEFAULT_MAX_ATTEMPTS):
        if half_order not in HALF_ORDERS:raise ValueError('Use X_FIRST or Z_FIRST')
        if type(max_attempts) is not int or max_attempts<1:raise ValueError('Positive integer max_attempts required in B2')
        if getattr(noise,'transport_events',()):raise ValueError('Transport remains empty')
        self.noise=noise;self.params=getattr(noise,'params',None);self.timed=isinstance(noise,IdleIonNoise)
        if timing_policy not in ('SEQ','PAR_MR') or (not self.timed and timing_policy!='SEQ'):raise ValueError('Unsupported timing policy')
        self.policy=timing_policy if self.timed else 'timing-free'
        self.half_order=half_order;self.order=HALF_ORDERS[half_order];self.max_attempts=max_attempts
        self.segments={};self.locations=[]
        for hi,kind in enumerate(self.order):
            basis='+' if kind=='X' else '0';base=preparation_ops(basis)
            if self.timed:
                ops,meta,_,duration,readout=insert_preparation_idles(base,self.params,timing_policy,live_on_entry=DATA)
            else:ops=base;meta=[{'part':'initial_resets' if i<8 else 'encoder' if i<16 else 'verification' if i<19 else 'verification_readout','original_op_index':i} for i in range(len(base))];duration=readout=None
            for attempt in range(max_attempts):self._segment(('P',hi,attempt),ops,meta,duration,readout,basis)
            base=extraction_ops(kind)
            if self.timed:ops,meta,duration,readout=extraction_idles(base,self.params,timing_policy)
            else:ops=base;meta=[{'part':'coupling' if i<7 else 'ancilla_readout','original_op_index':i} for i in range(len(base))];duration=readout=None
            self._segment(('C',hi,0),ops,meta,duration,readout,basis)

    def _segment(self,key,ops,meta,duration,readout,basis):
        self.segments[key]={'ops':ops,'metadata':meta,'duration':duration,'readout':readout,'basis':basis}
        self.locations.extend((key,i,op) for i,op in enumerate(ops))

    def channel(self,location):
        key,i,op=location
        if op[0]=='IDLE':
            p=idle_probability(op[2],self.params.T2)
            return {'Z':p} if p else {}
        return self.noise.channel((key[0],key[1],i,op))

    def run(self,*,faults=None,incoming=(0,0),basis='0',tableau=False,force_rejections=(0,0)):
        """force_rejections is test-only branch conditioning, not physical noise.

        Corrections are ideal Pauli-frame updates, once per extracted component.
        Aborts stop the round and are logical failures, with costs already spent.
        """
        faults=faults or {};x,z=incoming
        sim=None
        if tableau:
            sim=stim.TableauSimulator();ft.prepare_logical(sim,basis)
            for q in DATA:
                if x>>q&1:sim.x(q)
                if z>>q&1:sim.z(q)
        log={'aborted':False,'verification':[],'syndromes':{},'corrections':{},'executed':[]}
        def run(key):
            nonlocal x,z,sim
            log['executed'].append(key);ops=self.segments[key]['ops'];f=faults.get(key,{})
            if tableau:sim,bits=run_tableau(ops,f,sim=sim)
            else:x,z,bits=scalar.run_ops_frame(x,z,ops,f)
            return bits
        for hi,kind in enumerate(self.order):
            accepted=False
            for a in range(self.max_attempts):
                bit=run(('P',hi,a))['verify'];reject=bool(bit or a<force_rejections[hi])
                log['verification'].append({'half':kind,'attempt':a+1,'outcome':bit,'rejected':reject})
                if not reject:accepted=True;break
            if not accepted:log['aborted']=True;break
            bits=run(('C',hi,0));s=syndrome(bits);fix=ft.standard_decode(s)
            log['syndromes'][kind]=s;log['corrections'][kind]=list(fix)
            for q in fix:
                if tableau:
                    if kind=='X':sim.x(q)
                    else:sim.z(q)
                elif kind=='X':x^=1<<q
                else:z^=1<<q
        log['residual_signature']=tableau_signature(sim,basis) if tableau else frame_signature(x,z,basis)
        if tableau:
            failed=not ft.output_is_correctable(sim.copy(),basis)
        else:failed=bool(scalar.logical_failure(x,z)[0 if basis=='0' else 1])
        log['logical_failure']=bool(log['aborted'] or failed)
        return (sim if tableau else (x&MASK,z&MASK)),log

    def resources(self,executed):
        counts=Counter();seconds=0.;exposure={str(q):{'seconds':0.,'summed_Z_probability':0.} for q in DATA}
        for key in executed:
            segment=self.segments[key]
            seconds+=segment['duration'] or 0.
            for op in segment['ops']:
                kind=op[0];typ='reset' if kind in ('R','RX') else '2q' if kind=='CX' else 'measurement' if kind in ('M','MX') else 'idle'
                counts[typ]+=1
                if kind=='IDLE' and op[1] in DATA:
                    e=exposure[str(op[1])];e['seconds']+=op[2];e['summed_Z_probability']+=idle_probability(op[2],self.params.T2)
        return {'location_counts':dict(counts),'duration_seconds':seconds if self.timed else None,'ion_count':15,'idle_exposure_per_data_qubit':exposure if self.timed else None}
