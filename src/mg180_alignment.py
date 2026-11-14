"""Bounded independent secondary motion with mass-preserving finite geometry beam."""
import numpy as np
from scipy.special import logsumexp
from mg172_corrected import Corrected
CLASSES=['absent','single','Xonly','Yonly','multiple','displaced']
MOTIONS=['S','M'];PAIRS=[(k,m) for k in CLASSES for m in MOTIONS]
class AlignedTracker(Corrected):
    def __init__(self,*args,secondary_keep=8):
        super().__init__(*args);self.secondary_keep=secondary_keep;self.offsets=np.r_[np.arange(-24,0),np.arange(1,25)];self.displaced_norm=np.array([[self.inner(self.mean('displaced',0,si,(int(delta),0)),self.mean('displaced',0,si,(int(delta),0))) for delta in self.offsets] for si in range(13)])
        root=dict(kind='displaced',center=self.seed,si=4,second=None,mean=self.mean('single',self.seed,4,None),log=-np.log(6));self.beam.append(root);self.beam=[dict(p,motion=m,log=-np.log(12)) for p in self.beam for m in MOTIONS]
    def mean(self,kind,k,si,second):
        if kind!='displaced':return super().mean(kind,k,si,second)
        f=self.ft[si]*self.phase[k];f=f.copy()
        if second is not None:f[1]*=np.exp(-2j*np.pi*np.arange(251)*second[0]/500)
        return f
    def offset_prior(self,parent):
        p=np.full(48,1/48)
        if parent['kind']=='displaced' and parent['second'] is not None:
            old=parent['second'][0];ix=np.where(abs(self.offsets-old)<=1)[0];p[:]=0;p[ix]=.5/len(ix);p[np.where(self.offsets==old)[0][0]]+=.5
        return np.where(p>0,np.log(np.maximum(p,1e-300)),-np.inf)

    def step(self,field):
        y=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);rho=0 if self.previous is None else self.rho;innovation=y if self.previous is None else y-rho*self.previous;variance=1 if self.previous is None else 1-rho*rho;old=self.beam;resets=0;ends={p:{} for p in PAIRS};merge=0;shiftcache={};eligiblecache={};parentinfo=[];cache={}
        for parent in old:
            keymean=None if rho==0 else(parent['kind'],parent['center'],parent['si'],str(parent['second']))
            if keymean in cache:baseq,corr=cache[keymean]
            else:
                v=innovation+rho*parent['mean'];baseq=float(self.inner(v,v)/variance);corr={k:np.fft.irfft((v[None,:,:]*np.conj(K)).sum(axis=1),n=500,axis=1) for k,K in self.kernels.items()};cache[keymean]=(baseq,corr)
            parentinfo.append((baseq,corr))
            for motion in MOTIONS:
                logmode=np.log(.98 if motion==parent['motion'] else .02);pair=('absent',motion);key=(parent['center'],parent['si']);score=parent['log']+logmode+np.log(.98 if parent['kind']=='absent' else .004)-.5*baseq
                if key in ends[pair]:merge+=1;ends[pair][key]=float(np.logaddexp(ends[pair][key],score))
                else:ends[pair][key]=score
                locations=[key] if motion=='S' else [((parent['center']+d)%500,si) for d in range(-2,3) for si in range(max(0,parent['si']-1),min(13,parent['si']+2))]
                for k,si in locations:
                    key=(k,si);base=-.5*baseq+corr['single'][si,k]/variance-.5*self.norms[si]/variance
                    for kind in('single','Xonly','Yonly'):
                        pair=(kind,motion);score=parent['log']+logmode+np.log(.98 if kind==parent['kind'] else .004)-np.log(len(locations))-.5*baseq+corr[kind][si,k]/variance-.5*self.kindnorm[kind][si]/variance
                        if key in ends[pair]:merge+=1;ends[pair][key]=float(np.logaddexp(ends[pair][key],score))
                        else:ends[pair][key]=score
                    if k not in eligiblecache:
                        distance=np.minimum((np.arange(500)-k)%500,(k-np.arange(500))%500)*.2;eligiblecache[k]=np.broadcast_to(distance[None,:]>self.model.x.background_exclusion,(13,500))
                    eligible=eligiblecache[k];N=int(eligible.sum())
                    if parent['kind']=='multiple' and parent['second'] is not None:
                        oldk,oldsi=parent['second'];choices=[(ss,kk%500) for ss in range(max(0,oldsi-1),min(13,oldsi+2)) for kk in range(oldk-2,oldk+3) if eligible[ss,kk%500]]
                        priorgrid=np.full((13,500),-np.inf)
                        if choices:
                            probs={q:.5/len(choices) for q in choices}
                            if eligible[oldsi,oldk]:probs[(oldsi,oldk)]=probs.get((oldsi,oldk),0)+.5
                            normal=sum(probs.values())
                            for q,v in probs.items():priorgrid[q]=np.log(v/normal)
                        else:
                            priorgrid[eligible]=-np.log(N);resets+=1
                    else:
                        priorgrid=np.where(eligible,-np.log(N),-np.inf)

                    if key not in shiftcache:shiftcache[key]=np.roll(self.cross[si],k,axis=1)
                    delta=(corr['single']-shiftcache[key]-.5*self.norms[:,None])/variance;grid=parent['log']+logmode+np.log(.98 if parent['kind']=='multiple' else .004)-np.log(len(locations))+priorgrid+base+delta;grid=np.where(eligible,grid,-np.inf);pair=('multiple',motion)
                    if key in ends[pair]:merge+=1;ends[pair][key]=np.logaddexp(ends[pair][key],grid)
                    else:ends[pair][key]=grid
                    lp=self.offset_prior(parent);score=parent['log']+logmode+np.log(.98 if parent['kind']=='displaced' else .004)-np.log(len(locations))+lp-.5*baseq+(corr['Xonly'][si,k]+corr['Yonly'][si,(k+self.offsets)%500]-.5*self.displaced_norm[si])/variance;pair=('displaced',motion)
                    if key in ends[pair]:merge+=1;ends[pair][key]=np.logaddexp(ends[pair][key],score)
                    else:ends[pair][key]=score
        pairlogs=[];chosen={} 
        for pair in PAIRS:
            kind,motion=pair
            if kind=='multiple':
                masses=[];options=[]
                for (k,si),grid in ends[pair].items():
                    masses.append(float(logsumexp(grid)));flat=grid.ravel();indices=np.argpartition(flat,-self.secondary_keep)[-self.secondary_keep:]
                    for idx in indices:
                        ss,kk=np.unravel_index(int(idx),grid.shape);options.append(dict(kind=kind,motion=motion,center=k,si=si,second=(int(kk),int(ss)),merged_log=float(flat[idx])))
                pairlogs.append(float(logsumexp(masses)));chosen[pair]=sorted(options,key=lambda a:a['merged_log'],reverse=True)[:self.secondary_keep]
            elif kind=='displaced':
                masses=[];options=[]
                for (k,si),grid in ends[pair].items():
                    masses.append(float(logsumexp(grid)));indices=np.argpartition(grid,-self.secondary_keep)[-self.secondary_keep:]
                    for idx in indices:options.append(dict(kind=kind,motion=motion,center=k,si=si,second=(int(self.offsets[idx]),0),merged_log=float(grid[idx])))
                pairlogs.append(float(logsumexp(masses)));chosen[pair]=sorted(options,key=lambda a:a['merged_log'],reverse=True)[:self.secondary_keep]
            else:
                options=[dict(kind=kind,motion=motion,center=k,si=si,second=None,merged_log=float(value)) for (k,si),value in ends[pair].items()];pairlogs.append(float(logsumexp([a['merged_log'] for a in options])));chosen[pair]=sorted(options,key=lambda a:a['merged_log'],reverse=True)[:(1 if kind=='absent' else 2)]
        Z=float(logsumexp(pairlogs));classlogs=[float(logsumexp(pairlogs[2*i:2*i+2])) for i in range(len(CLASSES))];weights={k:float(np.exp(v-Z)) for k,v in zip(CLASSES,classlogs)};pairweights={k+'/'+m:float(np.exp(v-Z)) for (k,m),v in zip(PAIRS,pairlogs)};beam=[];fractions={}
        for pair,full in zip(PAIRS,pairlogs):
            selected=float(logsumexp([a['merged_log'] for a in chosen[pair]]));fractions[pair[0]+'/'+pair[1]]=float(np.exp(selected-full))
            for a in chosen[pair]:
                state=dict(a);state['mean']=np.zeros((2,251),complex) if a['kind']=='absent' else self.mean(a['kind'],a['center'],a['si'],a['second']);state['log']=a['merged_log']+full-selected-Z;beam.append(state)
        for pair in PAIRS:np.testing.assert_allclose(sum(np.exp(a['log']) for a in beam if(a['kind'],a['motion'])==pair),pairweights[pair[0]+'/'+pair[1]],rtol=1e-9,atol=1e-12)
        winner=max(weights,key=weights.get);best=max([a for pair in PAIRS if pair[0]==winner for a in chosen[pair]],key=lambda a:a['merged_log']);mu=np.zeros((2,251),complex) if winner=='absent' else self.mean(winner,best['center'],best['si'],best['second']);incoming=[]
        for parent,(baseq,corr) in zip(old,parentinfo):
            allowed=(best['center']==parent['center'] and best['si']==parent['si']) if winner=='absent' or best['motion']=='S' else min((best['center']-parent['center'])%500,(parent['center']-best['center'])%500)<=2 and abs(best['si']-parent['si'])<=1
            if not allowed:continue
            err=innovation-(mu-rho*parent['mean']);q=float(self.inner(err,err)/variance);score=parent['log']+np.log(.98 if winner==parent['kind'] else .004)+np.log(.98 if best['motion']==parent['motion'] else .02)-.5*q
            if winner!='absent' and best['motion']=='M':score-=np.log(5*(min(13,parent['si']+2)-max(0,parent['si']-1)))
            if winner=='multiple':
                eligible=eligiblecache[best['center']]
                if parent['kind']=='multiple' and parent['second'] is not None:
                    oldk,oldsi=parent['second'];choices=[(ss,kk%500) for ss in range(max(0,oldsi-1),min(13,oldsi+2)) for kk in range(oldk-2,oldk+3) if eligible[ss,kk%500]]
                    if choices:
                        probs={q:.5/len(choices) for q in choices}
                        if eligible[oldsi,oldk]:probs[(oldsi,oldk)]=probs.get((oldsi,oldk),0)+.5
                        ss,kk=best['second'][1],best['second'][0]
                        if (ss,kk) not in probs:continue
                        score+=np.log(probs[(ss,kk)]/sum(probs.values()))
                    else:score-=np.log(int(eligible.sum()))
                else:score-=np.log(int(eligible.sum()))
            if winner=='displaced':
                lp=self.offset_prior(parent)[np.where(self.offsets==best['second'][0])[0][0]]
                if not np.isfinite(lp):continue
                score+=lp
            incoming.append((score,q))
        qbest=max(incoming,key=lambda a:a[0])[1];one=[a for pair in PAIRS if pair[0]=='single' for a in chosen[pair]];spread=max(min((a['center']-b['center'])%500,(b['center']-a['center'])%500) for a in one for b in one);ws=max(abs(self.widths[a['si']]-self.widths[b['si']]) for a in one for b in one);label=winner if weights[winner]>=.95 and qbest<=self.qcut else 'uncertain'
        if label=='single' and(spread>4 or ws>.25):label='uncertain'
        singleZ=float(logsumexp(pairlogs[2:4]));singlemotion={m:float(np.exp(pairlogs[2+i]-singleZ)) for i,m in enumerate(MOTIONS)};self.beam=beam;self.previous=y;self.audit_snapshot=dict(old=old,innovation=innovation,variance=variance,ends=ends,parentinfo=parentinfo,pairlogs=pairlogs)
        return dict(label=label,winner=winner,center=best['center'],width=float(self.widths[best['si']]),representative_motion=best['motion'],class_weights=weights,pair_weights=pairweights,single_motion_weights=singlemotion,full_pair_logmass=pairlogs,full_class_logmass=classlogs,retained_pair_geometry_fraction=fractions,history_merge_operations=merge,secondary_kernel_resets=resets,selected_component_offset=(best['second'][0] if winner=='displaced' else None),innovation_q=qbest,qcut=self.qcut,single_center_spread=spread,single_width_spread=float(ws),beam=[dict(kind=a['kind'],motion=a['motion'],center=a['center'],width=float(self.widths[a['si']]),second=a['second'],merged_log=a['merged_log'],log_weight=a['log']) for a in beam])
