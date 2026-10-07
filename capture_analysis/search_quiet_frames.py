"""Bounded Eq.25 search of all pre-disturbance frames, with quiet row mean.
No frame stitching, decoding, whitening, calibrated SNR or detection claim.
Every overlapping eight-chirp window in frames0–33 and the quiet portion
of frame34 is searched on the finite3-D bank, followed by local refinement.
"""
from pathlib import Path
import json,hashlib,time
import h5py
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from analyze_captures import DATA,C,echo,read_profile
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'output';NC=8

class QuietBank:
    def __init__(self,p):
        self.p=p;self.u=np.tile(p['T_adc']+np.arange(p['samples'])/p['fs'],NC)
        self.t=(np.arange(NC)[:,None]*p['period']+(p['T_adc']+np.arange(p['samples'])/p['fs'])[None,:]).ravel()
        self.grids=[np.linspace(-1.2,1.2,49),np.linspace(.95*p['speed'],1.05*p['speed'],17),np.linspace(.05,1,20)]
        self.theta=np.stack(np.meshgrid(*self.grids,indexing='ij'),-1).reshape(-1,3)
        self.bounds=np.array([[g[0],g[-1]] for g in self.grids]);self.steps=np.array([g[1]-g[0] for g in self.grids])
        self.q=np.empty((len(self.theta),len(self.t)),np.complex64)
        for a in range(0,len(self.theta),256):self.q[a:a+256]=self.template(self.theta[a:a+256])
        print('Bank:',self.q.shape,flush=True)
    def template(self,theta):
        theta=np.atleast_2d(theta);x0,v,y=[theta[:,i,None] for i in range(3)]
        x=v*self.t-x0;rr=x*x+y*y
        d=rr/(np.sqrt((x*v)**2+(C*C-v*v)*rr)+x*v);tau=2*d;xs=x-v*d;r=np.hypot(xs,y)
        vr=v*xs/r;td=2*vr/(C+vr);s=self.p['freq_slope'];f0=self.p['f_start']
        f=-s*tau-(f0+s*(self.u-tau))*td
        phi=2*np.pi*(-f0*tau-s*self.u*tau+.5*s*tau*tau)
        hp=(1j*f/(175e3+1j*f))*(1j*f/(350e3+1j*f))
        # Full signed IQ hypotheses, not an assumed one-sided ADC visibility gate.
        valid=(abs(f)<self.p['if_max'])&(self.u>=tau)
        q=np.exp(1j*phi)*hp/r**2*valid
        # Record concentration diagnostically; do not exclude a physically
        # possible close passage merely because only a few chirps dominate.
        powers=np.sum(abs(q.reshape(-1,NC,self.p['samples']))**2,axis=2)
        eff=np.divide(powers.sum(1)**2,np.sum(powers**2,axis=1),out=np.zeros(len(q)),where=np.sum(powers**2,axis=1)>0)
        eligible=(valid.mean(1)>=.8)
        norm=np.linalg.norm(q,axis=1);return np.divide(q,norm[:,None],out=np.zeros_like(q),where=(eligible&(norm>0))[:,None])
    def coarse(self,z,starts,variance):
        # Save best template for each tested window, both stored-IQ orientations.
        best=np.zeros(len(starts));indices=np.zeros(len(starts),int);signs=np.zeros(len(starts),int)
        for a in range(0,len(starts),128):
            selected=np.stack([z[s:s+NC].ravel() for s in starts[a:a+128]]).T.astype(np.complex64)
            for sign in [0,1]:
                inner=(self.q.conj() if sign==0 else self.q)@selected
                power=abs(inner)**2/variance
                ix=power.argmax(0);val=power[ix,np.arange(len(ix))]
                use=val>best[a:a+len(ix)]
                best[a:a+len(ix)][use]=val[use];indices[a:a+len(ix)][use]=ix[use];signs[a:a+len(ix)][use]=sign
        return best,indices,signs
    def refine(self,data,seed,sign,variance):
        origin=self.bounds.mean(1);scale=np.array([.02,30.,.02]);bounds=(self.bounds-origin[:,None])/scale[:,None]
        def objective(w):
            q=self.template(origin+w*scale)[0]
            if sign:q=q.conj()
            return -abs(np.vdot(q,data))**2/variance
        result=minimize(objective,(seed-origin)/scale,method='Nelder-Mead',bounds=bounds,options={'maxiter':600,'xatol':1e-6,'fatol':1e-6})
        theta=origin+scale*result.x;q=self.template(theta)[0]
        if sign:q=q.conj()
        alpha=np.vdot(q,data);return float(-result.fun),theta,q*alpha,bool(result.success)

def selected_starts(z):
    energy=np.sum(abs(z)**2,axis=1)
    return np.arange(len(z)-NC+1),energy

def main():
    OUT.mkdir(exist_ok=True);starttime=time.monotonic()
    sourcepath=next(DATA.glob('2025*6607*.nc'));p=read_profile(sourcepath)
    with h5py.File(OUT/'quiet_subtracted_6607.h5') as sub:
        background=sub['background_complex_counts'][:];cutoff=int(sub.attrs['quiet_chirp_stop_exclusive'])
    with threadpool_limits(limits=4),h5py.File(sourcepath,'r') as raw,h5py.File(OUT/'quiet_matched_filter_6607.h5','w') as out:
        out.attrs['generator']='capture_analysis/search_quiet_frames.py';out.attrs['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        out.attrs['source_file']=sourcepath.name;out.attrs['source_id']=raw.attrs['id'];out.attrs['profile_json']=json.dumps(p)
        out.attrs['status']='Bounded measured search; candidate only, not a confirmed pellet echo'
        out.attrs['scope']='Every overlapping8-chirpwindow in frames0–33 and frame34chirps0–1628. No cross-frame windows. Two IQ signs. No whitening. Finite3Dbank,20localrefinements.'
        out.attrs['coordinate_definition']='x0 is radar longitudinal coordinate relative to projectile x=0 at first chirp of EACH tested window. y0 transverse distance; v0 positive speed. Unknown launch time absorbed into local x0. Not a global surveyed coordinate.'
        out.attrs['template_assumptions']='Eq25 phase,1/R2 amplitude,min175/350kHzHPFs,abs(IF)<10MHz, both IQ orientations, >=80%validsamples. No arbitrary1mtrajectorygate. ADC sideband unknown. Unknown complex amplitude; PEC frequency dependence omitted in normalized search.'
        out.create_dataset('quiet_background_complex_counts',data=background)
        # Last entirely quiet frame provides independent variance estimate;
        # score normalized by scalar sample variance, colored clutter not whitened.
        train=raw['radar_cube'][0,33,:,0,:].astype(np.complex128)-background
        variance=float(np.mean(abs(train)**2));out.attrs['reference_variance_counts2']=variance
        bank=QuietBank(p)
        for name,g in zip(['x0_grid_m','v0_grid_m_s','y0_grid_m'],bank.grids):out.create_dataset(name,data=g)
        reports=[];top=[]
        for frame in range(35):
            original=raw['radar_cube'][0,frame,:,0,:];z=original.astype(np.complex128)-background
            if frame==34:z=z[:cutoff]
            starts,energy=selected_starts(z)
            scores,indices,signs=bank.coarse(z,starts,variance)
            group=out.create_group(f'frame_{frame}');group.attrs['source_rx0_sha256']=hashlib.sha256(original.tobytes()).hexdigest()
            for name,val in [('window_starts',starts),('coarse_score',scores),('best_theta_index',indices),('iq_conjugated',signs),('residual_chirp_energy_counts2',energy)]:group.create_dataset(name,data=val)
            order=np.argsort(scores)[::-1];chosen=[]
            for i in order:
                if all(abs(int(starts[i])-int(starts[j]))>=NC for j in chosen):chosen.append(i)
                if len(chosen)==3:break
            for i in chosen:
                s=int(starts[i]);top.append((float(scores[i]),frame,s,bank.theta[indices[i]],int(signs[i]),z[s:s+NC].copy()))
            report=dict(frame=frame,quiet_chirps=len(z),tested_windows=len(starts),maximum_coarse_score=float(scores.max()),start_chirp=int(starts[order[0]]))
            reports.append(report);print(json.dumps(report),flush=True);out.flush()
        refined=[]
        for score,frame,s,seed,sign,z in sorted(top,key=lambda r:r[0],reverse=True)[:20]:
            val,theta,fitted,success=bank.refine(z.ravel(),seed,sign,variance)
            report=dict(frame=frame,start_chirp=s,chirps=NC,score=val,coarse_score=score,theta=theta.tolist(),iq_conjugated=bool(sign),optimizer_success=success,explained_residual_energy_fraction=val*variance/np.sum(abs(z)**2))
            refined.append((report,z,fitted.reshape(NC,p['samples'])));print('REFINED',json.dumps(report),flush=True)
        refined.sort(key=lambda r:r[0]['score'],reverse=True)
        out.attrs['frame_reports_json']=json.dumps(reports)
        for i,(report,z,fit) in enumerate(refined):
            gr=out.create_group(f'candidate_{i}');gr.attrs['report_json']=json.dumps(report);gr.create_dataset('residual_complex_counts',data=z);gr.create_dataset('fitted_complex_counts',data=fit)
        # Injection check in a real quiet window, same bank/refinement pipeline.
        reference=raw['radar_cube'][0,34,1500:1508,0,:].astype(np.complex128)-background
        theta=bank.theta[np.flatnonzero(np.linalg.norm(bank.q,axis=1)>.9)[len(bank.q)//5 % np.count_nonzero(np.linalg.norm(bank.q,axis=1)>.9)]]
        q=bank.template(theta)[0];injected=reference+np.sqrt(1000*variance)*q.reshape(NC,p['samples'])
        sc,ix,sg=bank.coarse(injected,np.array([0]),variance)
        val,rec,fit,success=bank.refine(injected.ravel(),bank.theta[ix[0]],int(sg[0]),variance)
        recovered=bank.template(rec)[0];recovered=recovered.conj() if sg[0] else recovered
        coherence=float(abs(np.vdot(recovered,q))**2)
        out.attrs['injection_check_json']=json.dumps(dict(injected_score=1000,recovered_score=val,truth=theta.tolist(),recovered=rec.tolist(),coherence=coherence))
        assert coherence>.9,'Injected template not recovered'
        out.attrs['elapsed_seconds']=time.monotonic()-starttime
        print('INJECTION',out.attrs['injection_check_json'],flush=True)
        print('BEST',json.dumps(refined[0][0]),flush=True)
        print('Written',out.filename,'elapsed',out.attrs['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
