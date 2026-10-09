import assert from 'node:assert/strict';
import fs from 'node:fs';
import {Heatmap} from '../web/lab/plots.mjs';
class Element {
  constructor(){this.children=[];this.listeners=new Map();this.dataset={};this.style={};this.offsetWidth=30;this.classList={add(){},remove(){}};}
  append(...children){this.children.push(...children);}
  setAttribute(){}
  addEventListener(name,fn){this.listeners.set(name,fn);}
  setPointerCapture(id){this.capture=id;}
  hasPointerCapture(id){return this.capture===id;}
  releasePointerCapture(){this.capture=null;}
  querySelector(selector){return this.children.find(c=>c.className.split(' ').includes(selector.slice(1)));}
  querySelectorAll(){return this.children;}
}
globalThis.document={createElement:()=>new Element(),documentElement:new Element()};
globalThis.window={getSelection:()=>({removeAllRanges(){}})};
const map=Object.create(Heatmap.prototype);
map.root=new Element();map.root.clientWidth=1176;
map.pixels={getBoundingClientRect:()=>({left:0,width:1000,height:200})};
map.config={x0:0,x1:100};
const event=(x,target)=>({button:0,pointerId:1,clientX:x,target,preventDefault(){},stopPropagation(){}});
function add(label,range,resizable=true){
  let commits=0;
  map.enableWindow({label,getRange:()=>range,enabled:()=>true,limit:()=>100,resizable,preview:(a,b)=>{range[0]=a;range[1]=b;},commit:()=>commits++});
  const band=map.windows.at(-1).band;
  const drag=(target,delta,cancel=false)=>{
    band.listeners.get('pointerdown')(event(100,target));
    band.listeners.get('pointermove')(event(100+delta*10,target));
    band.listeners.get(cancel?'pointercancel':'pointerup')(event(100+delta*10,target));
  };
  return {range,band,drag,commits:()=>commits};
}
for(const label of ['Background','Analyze interval']){
  const w=add(label,[10,30]);
  w.drag(w.band,5);assert.deepEqual(w.range,[15,35]);
  w.drag(w.band.querySelector('.window-label'),-3);assert.deepEqual(w.range,[12,32]);
  w.drag(w.band.querySelector('.start'),2);assert.deepEqual(w.range,[14,32]);
  w.drag(w.band.querySelector('.stop'),4);assert.deepEqual(w.range,[14,36]);
  w.drag(w.band,10,true);assert.deepEqual(w.range,[14,36]);assert.equal(w.commits(),4);
}
const coherent=add('Coherent integration',[18,26],false);
assert(!coherent.band.querySelector('.start'));coherent.drag(coherent.band,2);assert.deepEqual(coherent.range,[20,28]);
map.drawWindow();
for(const {band} of map.windows){assert.equal(band.style.top,'37px');assert.equal(band.style.height,'200px');}
assert(Number(coherent.band.style.zIndex)>Number(map.windows[0].band.style.zIndex),'Coherent window remains above wider intervals');
const css=fs.readFileSync('web/lab/style.css','utf8').match(/\.analysis-window\s*\{([^}]+)\}/)[1];
assert(!/\b(height|top|pointer-events)\s*:/.test(css),'Analyze window inherits full-height background drag geometry');
console.log('PASS: identical body/label dragging, both resize edges, cancellation, full-height layout and coherent-window access');
