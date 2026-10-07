"""Measured Eq25 sweep: -7 to +7 km/s,100m/s; all strict quiet windows.
Base conda. Positive-speed canonical bank exactly represents negative-speed
mirror(x0,v0)->(-x0,-v0). Unknown absolute direction cannot be identified here.
MPS acceleration uses real matrix products, verified against CPU complex BLAS.
"""
from pathlib import Path
import hashlib,json,time
import h5py,numpy as np
from threadpoolctl import threadpool_limits
import torch
from search_quiet_frames import QuietBank,NC,OUT
from analyze_captures import DATA,read_profile

class SignedBank(QuietBank):
    def __init__(self,p):
        self.p=p;self.u=np.tile(p['T_adc']+np.arange(128)/p['fs'],NC)
        self.t=(np.arange(NC)[:,None]*p['period']+(p['T_adc']+np.arange(128)/p['fs'])[None,:]).ravel()
        self.grids=[np.linspace(-1.2,1.2,49),np.arange(0,7000.1,100),np.linspace(.05,1,20)]
        self.signed_velocity=np.arange(-7000,7000.1,100)
        self.theta=np.stack(np.meshgrid(*self.grids,indexing='ij'),-1).reshape(-1,3)
        self.bounds=np.array([[g[0],g[-1]] for g in self.grids]);self.steps=np.array([g[1]-g[0] for g in self.grids])
        self.q=np.empty((len(self.theta),1024),np.complex64)
        for a in range(0,len(self.theta),256):self.q[a:a+256]=self.template(self.theta[a:a+256])
        self.device='mps' if torch.backends.mps.is_available() else 'cpu'
        self.qr=torch.tensor(self.q.real.copy(),device=self.device)
        self.qi=torch.tensor(self.q.imag.copy(),device=self.device)
        # Negative-speed mirror is an EXACT model symmetry, including light time.
        tests=np.array([[.3,1400,.2],[-.9,7000,.8],[.5,0,.4]])
        mirror=tests*np.array([-1,-1,1]);np.testing.assert_allclose(self.template(tests),self.template(mirror),rtol=1e-12,atol=1e-12)
        rng=np.random.default_rng(20261006);d=(rng.normal(size=(1024,7))+1j*rng.normal(size=(1024,7))).astype(np.complex64)
        gpu=self.accelerated(d)
        for sign in [0,1]:
            q=self.q.conj() if sign==0 else self.q
            ref=abs(q[:256]@d)**2
            np.testing.assert_allclose(gpu[sign][:256].cpu().numpy(),ref,rtol=3e-4,atol=1e-4)
        self.validation_max_abs_error=float(max(abs(gpu[sign][:256].cpu().numpy()-abs((self.q.conj() if sign==0 else self.q)[:256]@d)**2).max() for sign in [0,1]))
        del gpu
        print('Canonical bank',self.q.shape,'signed velocities',self.signed_velocity[0],self.signed_velocity[-1],'step100;device',self.device,'CPU/MPS max error',self.validation_max_abs_error,flush=True)
    def accelerated(self,d):
        dr=torch.tensor(d.real.copy(),device=self.device);di=torch.tensor(d.imag.copy(),device=self.device)
        a=self.qr@dr;b=self.qi@di;c=self.qr@di;e=self.qi@dr
        return ((a+b)**2+(c-e)**2,(a-b)**2+(c+e)**2)
    def coarse(self,z,starts,variance):
        best=np.zeros(len(starts));indices=np.zeros(len(starts),int);signs=np.zeros(len(starts),int)
        for a in range(0,len(starts),256):
            d=np.stack([z[s:s+NC].ravel() for s in starts[a:a+256]]).T.astype(np.complex64)
            powers=self.accelerated(d)
            for sign in [0,1]:
                val,ix=torch.max(powers[sign],dim=0)
                val=val.cpu().numpy()/variance;ix=ix.cpu().numpy();end=a+len(ix);use=val>best[a:end]
                best[a:end][use]=val[use];indices[a:end][use]=ix[use];signs[a:end][use]=sign
            del powers
        return best,indices,signs
    def surface(self,z,variance):
        powers=self.accelerated(z.ravel().astype(np.complex64)[:,None])
        positive=torch.maximum(powers[0],powers[1]).cpu().numpy().ravel().reshape(49,71,20)/variance
        return np.concatenate([positive[::-1,1:,:][:,::-1,:],positive],axis=1)

def main():
    began=time.monotonic();outpath=OUT/'quiet_signed_velocity_6607.h5';temp=outpath.with_suffix('.running.h5')
    with h5py.File(OUT/'quiet_subtracted_6607.h5') as sub:background=sub['background_complex_counts'][:]
    with h5py.File(OUT/'quiet_cutoff_audit_6607.h5') as audit:cutoff=int(audit.attrs['disturbance_first_chirp'])
    file=next(DATA.glob('2025*6607*.nc'));p=read_profile(file)
    with threadpool_limits(limits=8),torch.no_grad(),h5py.File(file) as raw,h5py.File(temp,'w') as out:
        out.attrs['status']='Running finite bank; no confirmed pellet echo'
        out.attrs['script']='capture_analysis/search_signed_velocities.py';out.attrs['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        out.attrs['source_file']=file.name;out.attrs['source_id']=raw.attrs['id'];out.attrs['profile_json']=json.dumps(p);out.attrs['strict_cutoff_chirp']=cutoff
        out.attrs['scope']='Every overlapping8chirpwindow in frames0–33and34before1624;no frame stitching;RX0;bothIQsigns; no whitening. Same quiet complex mean as displayed plot.'
        out.attrs['velocity_symmetry']='Negative bank computed exactly by(x0,v0)->(-x0,-v0); signed along-path direction not identifiable with free x0 and isotropic geometry. Signed velocity is NOT radial velocity.'
        out.attrs['template_assumptions']='Eq25 phase,1/R2 amplitude,minimumHPFs,absIF<10MHz,unknown complex amplitude, no arbitrary1mtrajectorygate; ADC sideband/geometry unverified.'
        out.create_dataset('quiet_background_complex_counts',data=background)
        train=raw['radar_cube'][0,33,:,0,:].astype(np.complex128)-background;variance=float(np.mean(abs(train)**2));out.attrs['reference_variance_counts2']=variance
        bank=SignedBank(p);out.attrs['accelerator']=bank.device;out.attrs['CPU_MPS_max_abs_test_error']=bank.validation_max_abs_error
        for key,val in [('x0_grid_m',bank.grids[0]),('v0_grid_m_s',bank.signed_velocity),('canonical_v0_grid_m_s',bank.grids[1]),('y0_grid_m',bank.grids[2])]:out.create_dataset(key,data=val)
        reports=[];top=[]
        # Report the event-adjacent quiet frame first, then complete all controls.
        for frame in [34]+list(range(34)):
            original=raw['radar_cube'][0,frame,:,0,:];z=original.astype(np.complex128)-background
            if frame==34:z=z[:cutoff]
            starts=np.arange(len(z)-NC+1);scores,indices,signs=bank.coarse(z,starts,variance)
            gr=out.create_group(f'frame_{frame}');gr.attrs['source_rx0_sha256']=hashlib.sha256(original.tobytes()).hexdigest()
            for k,val in [('window_starts',starts),('coarse_score',scores),('canonical_theta_index',indices),('iq_conjugated',signs)]:gr.create_dataset(k,data=val)
            order=np.argsort(scores)[::-1];chosen=[]
            for i in order:
                if all(abs(int(starts[i])-int(starts[j]))>=NC for j in chosen):chosen.append(i)
                if len(chosen)==3:break
            for i in chosen:top.append((float(scores[i]),frame,int(starts[i]),bank.theta[indices[i]].copy(),int(signs[i]),z[starts[i]:starts[i]+NC].copy()))
            report=dict(frame=frame,tested_windows=len(starts),maximum_score=float(scores.max()),start_chirp=int(starts[order[0]]),canonical_theta=bank.theta[indices[order[0]]].tolist())
            reports.append(report);print(json.dumps(report),flush=True);out.flush()
            if frame==34:
                start=int(starts[order[0]]);data=z[start:start+NC];gr.create_dataset('best_signed_coarse_surface',data=bank.surface(data,variance),compression='gzip')
                val,theta,fit,success=bank.refine(data.ravel(),bank.theta[indices[order[0]]],int(signs[order[0]]),variance)
                local=dict(frame=34,start_chirp=start,chirps=NC,score=val,coarse_score=float(scores.max()),theta=theta.tolist(),mirror_theta=(theta*np.array([-1,-1,1])).tolist(),iq_conjugated=bool(signs[order[0]]),optimizer_success=success,explained_residual_energy_fraction=float(val*variance/np.sum(abs(data)**2)))
                gr.attrs['best_refined_report_json']=json.dumps(local);gr.create_dataset('best_residual_complex_counts',data=data);gr.create_dataset('best_fitted_complex_counts',data=fit.reshape(NC,128));print('STRICT FRAME34',json.dumps(local),flush=True)
        refined=[]
        for val,frame,s,seed,sign,z in sorted(top,key=lambda r:r[0],reverse=True)[:20]:
            score,theta,fit,success=bank.refine(z.ravel(),seed,sign,variance)
            report=dict(frame=frame,start_chirp=s,chirps=NC,score=score,coarse_score=val,theta=theta.tolist(),mirror_theta=(theta*np.array([-1,-1,1])).tolist(),iq_conjugated=bool(sign),optimizer_success=success,explained_residual_energy_fraction=float(score*variance/np.sum(abs(z)**2)))
            refined.append((report,z,fit.reshape(NC,128)));print('REFINED',json.dumps(report),flush=True)
        refined.sort(key=lambda r:r[0]['score'],reverse=True)
        for i,(report,z,fit) in enumerate(refined):
            gr=out.create_group(f'candidate_{i}');gr.attrs['report_json']=json.dumps(report);gr.create_dataset('residual_complex_counts',data=z);gr.create_dataset('fitted_complex_counts',data=fit)
        # Broad-bank injected waveform check in measured quiet background.
        reference=raw['radar_cube'][0,34,1400:1408,0,:].astype(np.complex128)-background
        truth=np.array([.6,3500.,.35]);q=bank.template(truth)[0];injected=reference+np.sqrt(1000*variance)*q.reshape(NC,128)
        sc,ix,sg=bank.coarse(injected,np.array([0]),variance);score,rec,fit,success=bank.refine(injected.ravel(),bank.theta[ix[0]],int(sg[0]),variance)
        qr=bank.template(rec)[0];qr=qr.conj() if sg[0] else qr;coherence=float(abs(np.vdot(qr,q))**2)
        out.attrs['injection_check_json']=json.dumps(dict(truth=truth.tolist(),recovered=rec.tolist(),score=score,coherence=coherence));assert coherence>.9,'Broad-bank injection not recovered'
        out.attrs['frame_reports_json']=json.dumps(reports);out.attrs['elapsed_seconds']=time.monotonic()-began;out.attrs['status']='Complete finite-bank measured sweep; candidate only, not confirmed pellet echo'
        print('INJECTION',out.attrs['injection_check_json'],flush=True);print('BEST',json.dumps(refined[0][0]),flush=True);print('Elapsed',out.attrs['elapsed_seconds'],flush=True)
    temp.replace(outpath);print('Wrote',outpath,flush=True)
if __name__=='__main__':main()
