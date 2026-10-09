import fs from 'node:fs';import assert from 'node:assert/strict';
const p=JSON.parse(fs.readFileSync('web/lab/datasets/test63/metadata.json')).parameters;
const w=(await WebAssembly.instantiate(fs.readFileSync('web/lab/core.wasm'),{})).instance.exports;
const data=new Float32Array(2*225*125);data.fill(1);const ptr=w.allocate(data.length);new Float32Array(w.memory.buffer,ptr,data.length).set(data);
w.load(ptr,125,225,125,p.fs,p.T_adc,25.37e-6,p.f_start,p.freq_slope);w.release(ptr,data.length);w.background(0,2);
for(const [nc,nr,nv,na,groups] of [[1,329,53,2,2],[2,329,128,3,6],[4,329,255,10,30],[8,329,509,37,185]]){assert(w.prepare_power(4,nc,0,2,0)>0);assert.equal(w.radial_begin(.001,3,0,600,0,1e6,.05,32000000),groups);w.radial_info();const info=Array.from(new Float32Array(w.memory.buffer,w.result_ptr(),w.result_len()));assert.deepEqual(info.slice(0,3),[nr,na,nv]);}
for(const nc of [16,32]){w.prepare_power(4,nc,0,2,0);assert(w.radial_begin(.001,3,0,600,0,1e6,.05,32000000)<0);}
console.log('PASS: complexity table matches shipped Wasm planner; 16/32 full banks exceed cap');
