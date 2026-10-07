"""RX0 complex spectrum of EACH stored chirp, no range conversion or background subtraction.

Run with conda run -n base python capture_analysis/plot_chirp_spectra.py.
Hann window, 128 point complex FFT, signed frequency. All 16384 chirps per shot.
"""
from pathlib import Path
import hashlib,json
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_captures import DATA

ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'output'

def main():
    for shot,frames,name in [('6606',[43,44,45,46],'parallel'),('6607',[34,35,36,37],'perpendicular')]:
        with h5py.File(OUT/'raw_receiver_iq.h5') as prior:
            group=prior[name];p=json.loads(group.attrs['profile_json']);source_id=group.attrs['source_id']
            marked=int(group.attrs['frame']);onset=int(group.attrs['data_selected_onset_chirp'])
        with h5py.File(DATA/p['file']) as raw:
            blocks=[raw['radar_cube'][0,f,:,0,:] for f in frames]
        z=np.concatenate(blocks).astype(np.complex128)
        window=np.hanning(z.shape[1]);f=np.fft.fftshift(np.fft.fftfreq(z.shape[1],1/p['fs']))
        transform=np.fft.fftshift(np.fft.fft(z*window[None,:],axis=1),axes=1)/window.sum()
        power=abs(transform)**2;reference=float(power.max())
        db=10*np.log10(np.maximum(power,np.finfo(float).tiny)/reference)
        suffix=f'{shot}_frames_'+ '_'.join(map(str,frames))
        output=OUT/f'chirp_spectra_{suffix}.h5'
        with h5py.File(output,'w') as h:
            h.attrs['script']='capture_analysis/plot_chirp_spectra.py'
            h.attrs['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            h.attrs['source_file']=p['file'];h.attrs['source_id']=source_id
            h.attrs['frames']=frames;h.attrs['receiver']=0;h.attrs['fs_hz']=p['fs']
            h.attrs['chirp_period_s_assumed']=p['period'];h.attrs['frame_gaps']='omitted; secondary axis is concatenated chirp time, not elapsed time across unknown frame gaps'
            h.attrs['processing']='Each complex chirp separately; 128-sample Hann; FFT128; coherent-amplitude-normalized FFT; no subtraction, mean removal, conjugation, masking or whitening'
            h.attrs['reference_spectral_power_counts2']=reference
            h.attrs['source_rx0_sha256_json']=json.dumps({str(frame):hashlib.sha256(block.tobytes()).hexdigest() for frame,block in zip(frames,blocks)})
            h.create_dataset('complex_counts',data=z,compression='gzip')
            h.create_dataset('frequency_hz',data=f).attrs['units']='signed frequency of stored I+iQ, not inferred range or Doppler alone'
            h.create_dataset('spectral_power_counts2',data=power,compression='gzip')
            h.create_dataset('spectral_power_db_relative_to_peak',data=db.astype(np.float32),compression='gzip')
            h.create_dataset('chirp_index',data=np.arange(len(z)))
            h.create_dataset('concatenated_time_s',data=np.arange(len(z))*p['period'])
            h.create_dataset('source_frame_index',data=np.repeat(frames,len(blocks[0])))
        fig,ax=plt.subplots(figsize=(12,7),layout='constrained')
        freq_edges=(np.arange(len(f)+1)-len(f)/2-.5)*p['fs']/len(f)/1e6
        chirp_edges=np.arange(len(z)+1)-.5
        mesh=ax.pcolormesh(chirp_edges,freq_edges,db.T,cmap='magma',shading='flat',vmin=-70,vmax=0,rasterized=True)
        for i,frame in enumerate(frames):
            if i:ax.axvline(i*len(blocks[0])-.5,color='white',lw=.8,alpha=.8)
            ax.text((i+.5)*len(blocks[0]),.98,f'Frame {frame}',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=10,bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=2))
        # Preserve the same unlabeled visual boundary as the corresponding raw plot.
        ax.axvline(frames.index(marked)*len(blocks[0])+onset,color='white',ls='--',lw=.9,alpha=.8)
        ax.set(xlabel='Concatenated chirp index (frame gaps omitted)',ylabel='Signed sampled beat frequency (MHz)',xlim=(-.5,len(z)-.5),ylim=(-p['fs']/2e6,p['fs']/2e6))
        secondary=ax.secondary_xaxis('top',functions=(lambda x:x*p['period'],lambda t:t/p['period']))
        secondary.set_xlabel(f'Time since start of frame {frames[0]} (s; frame gaps omitted)')
        fig.suptitle(f'Shot {shot} — RX0, spectrum of each complex chirp\nFrames '+', '.join(map(str,frames))+f': all {len(z)} chirps; Hann window, {p["fs"]/len(f)/1e3:.1f} kHz bins',fontsize=14)
        bar=fig.colorbar(mesh,ax=ax,pad=.03,fraction=.05);bar.set_label('Spectral power (dB relative to maximum in these frames)')
        fig.supxlabel('Script: capture_analysis/plot_chirp_spectra.py · Raw complex I/Q; no background subtraction; frequency is not range',fontsize=7,color='gray')
        image=OUT/f'chirp_spectra_{suffix}.png';fig.savefig(image,dpi=200);plt.close(fig)
        print(shot,z.shape,'saved',image,'and',output,flush=True)
        assert power.shape==(16384,128) and np.isfinite(db).all()
        # Window-normalized Parseval, preserving the complex frequency sign.
        np.testing.assert_allclose(power.sum(1),len(f)*np.sum(abs(z*window)**2,axis=1)/window.sum()**2,rtol=1e-12)

if __name__=='__main__':main()
