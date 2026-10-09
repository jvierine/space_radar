// Verify the viewer's trajectory-completion event sequence without a browser.
import fs from 'node:fs';
import vm from 'node:vm';
import assert from 'node:assert/strict';
const elements=new Map(),listeners=new Map(),jobs=[];
function element(id){
 if(!elements.has(id))elements.set(id,{value:id==='#beamSteps'?'10':'',hidden:false,disabled:false,append(){},after(){},insertBefore(){},querySelector:element,querySelectorAll(){return []}});
 return elements.get(id);
}
const context={document:{createElement:()=>element(Symbol()),getElementById:element},window:{addEventListener:(name,fn)=>listeners.set(name,fn)},Heatmap:class{resetZoom(){}set(){}},Worker:class{postMessage(job){jobs.push(job)}},console};
const source=fs.readFileSync('web/lab/beam-ui.mjs','utf8').replace(/^import .*;\n/gm,'');
vm.runInNewContext(source,context);
const emit=(name,detail)=>listeners.get(name)({detail});
const fit=start=>({meta:{},match:{start,pulses:8}});
emit('fmcw-busy',true);
emit('fmcw-match',fit(750));
emit('fmcw-match',fit(875));
assert.equal(jobs.length,0,'Wait until the trajectory worker completes');
emit('fmcw-busy',false);
assert.equal(jobs.length,1);
assert.equal(jobs[0].match.start,875,'Use the final train from a comparison/scan');
assert.equal(jobs[0].steps,10);
emit('fmcw-busy',false);
assert.equal(jobs.length,1,'Do not repeat the automatic calculation');
emit('fmcw-invalidated');
emit('fmcw-busy',false);
assert.equal(jobs.length,1,'Do not calculate a stale fit');
emit('fmcw-match',fit(1000));
assert.equal(jobs.length,2,'A new completed fit calculates immediately');
console.log('PASS: RCS/beam calculation automatically runs once after the final fit; invalidated fits are excluded');
