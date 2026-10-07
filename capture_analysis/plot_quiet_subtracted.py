"""Subtract a pre-disturbance complex mean per fast-time bin; no decoding."""
from pathlib import Path
import hashlib, json
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_captures import DATA
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'output'

def main():
    script_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with h5py.File(OUT/'raw_receiver_iq.h5','r') as meta, h5py.File(OUT/'quiet_subtracted_6607.h5','w') as out:
        old=meta['perpendicular']; p=json.loads(old.attrs['profile_json'])
        frames=[34,35,36,37]; cutoff=int(old.attrs['data_selected_onset_chirp'])
        with h5py.File(DATA/p['file'],'r') as source:
            blocks=[source['radar_cube'][0,k,:,0,:] for k in frames]
        z=np.concatenate(blocks).astype(np.complex128)
        background=z[:cutoff].mean(axis=0)
        residual=z-background[None,:]
        out.attrs['script']='capture_analysis/plot_quiet_subtracted.py'
        out.attrs['script_sha256']=script_hash
        out.attrs['source_file']=old.attrs['source_file'];out.attrs['source_id']=old.attrs['source_id']
        out.attrs['frames']=frames;out.attrs['receiver']=0
        out.attrs['quiet_chirp_start_inclusive']=0;out.attrs['quiet_chirp_stop_exclusive']=cutoff
        out.attrs['processing']='Complex mean over quiet chirps, separately for each fast-time bin, subtracted from ALL chirps. No decoding, FFT, matched filter, whitening or power calibration.'
        out.attrs['frame_gap_note']='Frame gaps omitted; time assumes 13.71 us chirp spacing'
        out.attrs['source_rx0_sha256']=hashlib.sha256(np.concatenate(blocks).tobytes()).hexdigest()
        out.create_dataset('background_complex_counts',data=background)
        out.create_dataset('residual_complex_counts',data=residual,compression='gzip')
        out.create_dataset('chirp_index',data=np.arange(len(z)))
        out.create_dataset('concatenated_time_s',data=np.arange(len(z))*p['period'])
        out.create_dataset('fast_time_us',data=np.arange(z.shape[1])/p['fs']*1e6)
        print('quiet interval:',cutoff,'chirps; duration',cutoff*p['period'],'s')
        print('quiet residual RMS complex / real:',np.sqrt(np.mean(abs(residual[:cutoff])**2)),np.std(residual[:cutoff].real))
        for name,lo,hi,limit in [('full',0,len(z),32768),('transition',cutoff-250,cutoff+120,32768),('quiet_detail',0,cutoff+40,300)]:
            fig,ax=plt.subplots(figsize=(12,7),layout='constrained')
            m=ax.pcolormesh(np.arange(lo,hi+1)-.5,(np.arange(129)-.5)/p['fs']*1e6,residual[lo:hi].real.T,cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
            ax.axvline(cutoff,color='black',ls='--',lw=1)
            if name=='full':
                for i,f in enumerate(frames):
                    if i: ax.axvline(i*4096-.5,color='black',lw=.8)
                    ax.text((i+.5)*4096,.98,f'Frame {f}',transform=ax.get_xaxis_transform(),ha='center',va='top',bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=2))
            ax.set(xlabel='Concatenated chirp index (frame gaps omitted)',ylabel='Fast time within chirp (µs)',xlim=(lo-.5,hi-.5))
            a=ax.secondary_xaxis('top',functions=(lambda x:x*p['period'],lambda t:t/p['period']))
            a.set_xlabel('Time since start of frame 34 (s)')
            bar=fig.colorbar(m,ax=ax,pad=.03);bar.set_label('Real residual (ADC counts)')
            title='Shot 6607 — RX0, real part after quiet complex-background subtraction'
            subtitle=f'Background: mean of chirps 0–{cutoff-1}, per fast-time bin; no decoding'
            if name=='quiet_detail': subtitle+='\nColour clipped at ±300 counts to show weak changes'
            fig.suptitle(title+'\n'+subtitle,fontsize=13)
            fig.text(.01,.001,'Script: capture_analysis/plot_quiet_subtracted.py',fontsize=7,color='gray')
            fig.savefig(OUT/f'quiet_subtracted_6607_{name}.png',dpi=200)
            plt.close(fig)
            print(OUT/f'quiet_subtracted_6607_{name}.png')
if __name__=='__main__':main()
