"""Run Dawn/WebGPU numerical and fallback tests; save benchmark provenance in HDF5."""
from pathlib import Path
import hashlib,json,subprocess
import h5py
ROOT=Path(__file__).resolve().parents[1]
result=subprocess.run(['node','tools/validate_gpu_radial.mjs'],cwd=ROOT,stdout=subprocess.PIPE,text=True,check=True)
report=json.loads(result.stdout)
out=ROOT/'web/lab/qa/gpu_validation.h5'
out.parent.mkdir(parents=True,exist_ok=True)
with h5py.File(out,'w') as h:
    h.attrs['generator']='tools/run_gpu_validation.py + tools/validate_gpu_radial.mjs'
    h.attrs['gpu']=json.dumps(report['gpu'])
    for name in ['web/lab/core.wasm','web/lab/gpu-radial.mjs','web/lab/radial-backend.mjs']:
        h.attrs[name+'_sha256']=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    for row in report['records']:
        g=h.create_group(row['name'])
        for key,value in row.items():
            if key=='name':continue
            if isinstance(value,list):g.create_dataset(key,data=value)
            else:g.attrs[key]=value
print('PASS:',out)
