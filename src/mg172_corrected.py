"""Exact current secondary-prior sums within a class-mass-preserving geometry beam."""
import numpy as np
from scipy.special import logsumexp
from analyze_mg170 import Beam,CLASSES
class Corrected(Beam):
    def __init__(self,*args):
        super().__init__(*args);self.kernels={};self.kindnorm={};self.cross=[]
        for kind in('single','Xonly','Yonly'):
            f=self.ft.copy()
            if kind=='Xonly':f[:,1]=0
            if kind=='Yonly':f[:,0]=0
            self.kernels[kind]=np.einsum('fst,itf->isf',self.P,f);self.kindnorm[kind]=self.inner(f,f)
        for f in self.ft:self.cross.append(np.fft.irfft((f[None,:,:]*np.conj(self.kernels['single'])).sum(axis=1),n=500,axis=1))
        self.cross=np.array(self.cross)
        np.testing.assert_allclose([self.cross[i,i,0] for i in range(13)],self.norms,rtol=1e-10,atol=1e-8)
    def step(self,field):
        y=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);rho=0 if self.previous is None else self.rho;innovation=y if self.previous is None else y-rho*self.previous;variance=1 if self.previous is None else 1-rho*rho;old=self.beam;ends={k:{} for k in CLASSES};merge_count=0;shiftcache={};eligiblecache={};parentinfo=[]
        common=None
        for parent in old:
            v=innovation+rho*parent['mean']
            if rho==0 and common is not None:baseq,corr=common
            else:
                baseq=float(self.inner(v,v)/variance);corr={k:np.fft.irfft((v[None,:,:]*np.conj(K)).sum(axis=1),n=500,axis=1) for k,K in self.kernels.items()}
                if rho==0:common=(baseq,corr)
            locations=[((parent['center']+d)%500,si) for d in range(-2,3) for si in range(max(0,parent['si']-1),min(13,parent['si']+2))];parentinfo.append((baseq,corr,len(locations)))
            key=(parent['center'],parent['si']);value=parent['log']+np.log(.98 if parent['kind']=='absent' else .005)-.5*baseq
            if key in ends['absent']:merge_count+=1;ends['absent'][key]=float(np.logaddexp(ends['absent'][key],value))
            else:ends['absent'][key]=value
            for k,si in locations:
                key=(k,si);base=-.5*baseq+corr['single'][si,k]/variance-.5*self.norms[si]/variance
                for kind in('single','Xonly','Yonly'):
                    score=parent['log']+np.log(.98 if kind==parent['kind'] else .005)-np.log(len(locations))-.5*baseq+corr[kind][si,k]/variance-.5*self.kindnorm[kind][si]/variance
                    if key in ends[kind]:merge_count+=1;ends[kind][key]=float(np.logaddexp(ends[kind][key],score))
                    else:ends[kind][key]=score
                if k not in eligiblecache:
                    distance=np.minimum((np.arange(500)-k)%500,(k-np.arange(500))%500)*.2;eligiblecache[k]=np.broadcast_to(distance[None,:]>self.model.x.background_exclusion,(13,500))
                eligible=eligiblecache[k];N=int(eligible.sum())
                if key not in shiftcache:shiftcache[key]=np.roll(self.cross[si],k,axis=1)
                delta=(corr['single']-shiftcache[key]-.5*self.norms[:,None])/variance;loggrid=parent['log']+np.log(.98 if parent['kind']=='multiple' else .005)-np.log(len(locations))-np.log(N)+base+delta;loggrid=np.where(eligible,loggrid,-np.inf)
                if key in ends['multiple']:merge_count+=1;ends['multiple'][key]=np.logaddexp(ends['multiple'][key],loggrid)
                else:ends['multiple'][key]=loggrid
        classlogs=[];chosen={}
        for kind in CLASSES:
            if kind=='multiple':
                masses=[];options=[]
                for (k,si),grid in ends[kind].items():
                    masses.append(float(logsumexp(grid)));flat=grid.ravel();top=np.argpartition(flat,-2)[-2:]
                    for idx in top:
                        ss,kk=np.unravel_index(int(idx),grid.shape);options.append(dict(kind=kind,center=k,si=si,second=(int(kk),int(ss)),merged_log=float(flat[idx])))
                classlogs.append(float(logsumexp(masses)));chosen[kind]=sorted(options,key=lambda v:v['merged_log'],reverse=True)[:2]
            else:
                options=[dict(kind=kind,center=k,si=si,second=None,merged_log=float(value)) for (k,si),value in ends[kind].items()];classlogs.append(float(logsumexp([v['merged_log'] for v in options])));chosen[kind]=sorted(options,key=lambda v:v['merged_log'],reverse=True)[:(1 if kind=='absent' else 2)]
        total=float(logsumexp(classlogs));weights={k:float(np.exp(l-total)) for k,l in zip(CLASSES,classlogs)};beam=[];fractions={}
        for kind,full in zip(CLASSES,classlogs):
            selected=float(logsumexp([v['merged_log'] for v in chosen[kind]]));fractions[kind]=float(np.exp(selected-full))
            for state in chosen[kind]:
                st=dict(state);st['mean']=np.zeros((2,251),complex) if kind=='absent' else self.mean(kind,st['center'],st['si'],st['second']);st['log']=st['merged_log']+full-selected-total;beam.append(st)
        for kind in CLASSES:np.testing.assert_allclose(sum(np.exp(s['log']) for s in beam if s['kind']==kind),weights[kind],rtol=1e-9,atol=1e-12)
        winner=max(weights,key=weights.get);best=chosen[winner][0];mu=np.zeros((2,251),complex) if winner=='absent' else self.mean(winner,best['center'],best['si'],best['second']);incoming=[]
        for parent,(baseq,corr,nloc) in zip(old,parentinfo):
            allowed=(best['center']==parent['center'] and best['si']==parent['si']) if winner=='absent' else min((best['center']-parent['center'])%500,(parent['center']-best['center'])%500)<=2 and abs(best['si']-parent['si'])<=1
            if not allowed:continue
            q=float(self.inner(innovation-(mu-rho*parent['mean']),innovation-(mu-rho*parent['mean']))/variance);prior=parent['log']+np.log(.98 if winner==parent['kind'] else .005)
            if winner!='absent':prior-=np.log(nloc)
            if winner=='multiple':prior-=np.log(int(eligiblecache[best['center']].sum()))
            incoming.append((prior-.5*q,q))
        qbest=max(incoming,key=lambda v:v[0])[1];label=winner if weights[winner]>=.95 and qbest<=self.qcut else 'uncertain';one=chosen['single'];spread=max(min((a['center']-b['center'])%500,(b['center']-a['center'])%500) for a in one for b in one);ws=max(abs(self.widths[a['si']]-self.widths[b['si']]) for a in one for b in one)
        if label=='single' and(spread>4 or ws>.25):label='uncertain'
        self.beam=beam;self.previous=y
        result=dict(label=label,winner=winner,center=best['center'],width=float(self.widths[best['si']]),class_weights=weights,full_class_logmass=classlogs,retained_geometry_fraction=fractions,history_merge_operations=merge_count,innovation_q=qbest,qcut=self.qcut,single_center_spread=spread,single_width_spread=float(ws),beam=[dict(kind=s['kind'],center=s['center'],width=float(self.widths[s['si']]),second=s['second'],merged_log=s['merged_log'],log_weight=s['log']) for s in beam])
        self.audit_snapshot=dict(old=old,innovation=innovation,variance=variance,ends=ends,parentinfo=parentinfo,classlogs=classlogs)
        return result
