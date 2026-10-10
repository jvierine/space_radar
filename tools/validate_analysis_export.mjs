import assert from 'node:assert/strict';
import fs from 'node:fs';
import * as h5 from '../web/lab/vendor/h5wasm/hdf5_hl.js';
import {writeAnalysis} from '../web/lab/export-analysis.mjs';
const payload={url:'https://juha.no/fmcw/lab/',meta:{test:63},settings:{background:[0,125]},points:[{start:750,pulses:8,time:.019,kinematics:[.3,320,5e5],best:[.3,320,5e5,100],model:'radial',grid:{},beam:{rawNoisePower:[1,4,9,16],channelPower:[1,4,9,16].map(quiet=>({observed:100,quiet})),single:[10,20,30,40],peak:90,phases:[0,.1,.2,.3],projections:new Float32Array([1,2,3]),effectiveTime:.000144},rcs:Array.from({length:5},(_,i)=>({sigma:.001*(i+1),diameters:[3+i,4+i]}))}],result:{start:750,fit:new Float32Array([1,2,3,4]),observed:new Float32Array([2,3,4,5]),cube:new Float32Array([1,2,3]),best:[.3,320,5e5,100]}};
payload.result.rangeProfile=payload.points[0].rangeProfile={ranges:[.2,.3],observedPower:[4,9],noisePower:[.1,.1],score:[40,90],snrDb:[15.91,19.5],weights:[.5,0,.5,0,.5,0,.5,0],velocity:-320,acceleration:5e5,peakRange:.3,normalization:'fixed peak beam and motion'};
const bytes=await writeAnalysis(h5,payload);assert.deepEqual([...bytes.slice(0,8)],[137,72,68,70,13,10,26,10]);
fs.writeFileSync('/tmp/fmcw-export-validation.h5',bytes);
console.log('PASS: shipped browser HDF5 runtime writes analysis, all diameter roots, beam projections, grid and complex fit');

const {FS}=await h5.ready;FS.writeFile('/verify-profile.h5',bytes);
const file=new h5.File('/verify-profile.h5','r');
assert.deepEqual([...file.get('current_match/coherent_range_profile/ranges').value],[.2,.3]);
assert.deepEqual([...file.get('trains/chirp_750/coherent_range_profile/weights').value],payload.points[0].rangeProfile.weights);
file.close();FS.unlink('/verify-profile.h5');
console.log('PASS: current and historical coherent range profiles retained in HDF5');
