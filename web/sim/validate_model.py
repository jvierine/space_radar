"""Cross-language check against the memo's existing Eq.25 Python model.
Run with conda run -n base python web/sim/validate_model.py.
Reference samples and error metrics are saved in HDF5, not CSV.
"""
from pathlib import Path
import importlib.util,subprocess,json,hashlib
import numpy as np,h5py
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parents[1]/'capture_analysis/analyze_captures.py'
spec=importlib.util.spec_from_file_location('reference_echo',SOURCE);ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
js="""
import {PROFILES,simulate,state,timing,fft,spectrogram,alias} from './model.mjs';
const rows=[];
for(const profile of ['6606','6607'])for(const g of [{x0:.12,y0:.25,v:6285},{x0:-.35,y0:.15,v:-7000},{x0:5,y0:1,v:14000},{x0:1,y0:1,v:0}])for(const n of [1,10]){
  const p=PROFILES[profile],s=simulate(p,g,n,'receiver');
  if(s.points.length!==128*n)throw Error('sample count');
  for(let k=1;k<n;k++)if(Math.abs(s.records[k][0].t-s.records[k-1][0].t-timing(p).period)>1e-15)throw Error('chirp padding');
  const t=4e-6,u=4e-6,e=2e-9,derivative=(state(t+e,u+e,p,g).phase-state(t-e,u-e,p,g).phase)/(4*Math.PI*e);
  if(Math.abs(derivative-state(t,u,p,g).frequency)>10)throw Error('analytic frequency derivative');
  rows.push({profile,g,n,values:s.points.map(x=>[x.t,x.u,x.range,x.phase,x.frequency,x.i,x.q])});
}
for(const sign of [-1,1]) {
  const p=PROFILES['6607'],nfft=256,bin=16,f=sign*p.fs*bin/nfft;
  const record=Array.from({length:128},(_,j)=>({t:j/p.fs,k:0,i:Math.cos(2*Math.PI*f*j/p.fs),q:Math.sin(2*Math.PI*f*j/p.fs)}));
  const s=spectrogram({p,records:[record]});
  const peak=[...s.frames[0].power].indexOf(Math.max(...s.frames[0].power));
  if(peak!==nfft/2+sign*bin)throw Error('complex FFT frequency sign');
}
if(alias(7e6,12.5e6)!==-5.5e6 || alias(-8e6,12.5e6)!==4.5e6)throw Error('ADC alias');
console.log(JSON.stringify(rows));
"""
rows=json.loads(subprocess.run(['node','--input-type=module','-e',js],cwd=ROOT,check=True,capture_output=True,text=True).stdout)
(ROOT/'qa').mkdir(exist_ok=True)
out=ROOT/'qa/model_validation.h5';worst_phase=0.;worst_if=0.;worst_iq=0.
with h5py.File(out,'w') as h:
    h.attrs['reference_model']=str(SOURCE);h.attrs['reference_model_sha256']=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    h.attrs['web_model_sha256']=hashlib.sha256((ROOT/'model.mjs').read_bytes()).hexdigest()
    h.attrs['validator_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for idx,row in enumerate(rows):
        p=dict(f_start=78e9 if row['profile']=='6606' else 79.5e9,freq_slope=100e12 if row['profile']=='6606' else 60e12,fs=12.5e6,samples=128,T_adc=0,period=13.709999904632568e-6,if_max=10e6)
        g=row['g'];echo=ref.echo(p,[g['x0'],g['v'],g['y0']],row['n'],gate=False,phase_only=True)
        values=np.asarray(row['values']);iq=values[:,5]+1j*values[:,6];truth=echo['iq'].ravel()/echo['range'].ravel()**2
        expected=np.column_stack([echo[k].ravel() for k in ['t','u','range','phase','if_hz']])
        np.testing.assert_allclose(values[:,:5],expected,rtol=1e-12,atol=1e-8)
        np.testing.assert_allclose(iq,truth,rtol=1e-8,atol=1e-8)
        worst_phase=max(worst_phase,np.max(abs(values[:,3]-expected[:,3])))
        worst_if=max(worst_if,np.max(abs(values[:,4]-expected[:,4])))
        worst_iq=max(worst_iq,np.max(abs(iq-truth)))
        group=h.create_group(f'case_{idx:02d}');group.attrs['parameters_json']=json.dumps({k:v for k,v in row.items() if k!='values'})
        group.create_dataset('javascript_t_u_range_phase_if_i_q',data=values)
        group.create_dataset('python_reference_t_u_range_phase_if',data=expected)
        group.create_dataset('python_reference_iq',data=truth)
    h.attrs['status']='PASS: 16 parameter cases; memo phase/range/IF/receiver-IQ; padding; analytic derivative; complex FFT sign; sampling alias'
    h.attrs['max_phase_error_rad']=worst_phase;h.attrs['max_frequency_error_hz']=worst_if;h.attrs['max_relative_voltage_absolute_error']=worst_iq
print('PASS, 16 cases; phase error',worst_phase,'rad; IF error',worst_if,'Hz; complex voltage error',worst_iq)
print('Padding, instantaneous frequency derivative, complex FFT signs and ADC alias checks passed.')
print(out)
