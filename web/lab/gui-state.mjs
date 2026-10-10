const plots=new Map();
export function decodeState(url) {
  try {
    const s=JSON.parse(new URL(url).searchParams.get('gui'));
    return s?.version===1&&s.controls&&typeof s.controls==='object'?s:null;
  } catch {return null;}
}
export function encodeState(url,state) {
  const u=new URL(url);u.searchParams.set('gui',JSON.stringify({version:1,...state}));return u.href;
}
export const defaultState={"version":1,"controls":{"rx":"0","component":"0","rawScale":"2500","subScale":"100","span":"70","fftSub":"1","chirp":"3821","chirpSlider":"3821","traceSub":true,"algorithm":"fft","phaseLoss":"36.90426555198067","processingLossDb":"2","pulses":"8","computeBackend":"auto","receiver":false,"xMin":"0.001","xMax":"3","xN":"329","yMin":"0","yMax":"1000000","yN":"37","vMin":"0","vMax":"600","vN":"509","beamSteps":"10","rcsTemperature":"9000","rcsPower":"12","rcsTxGain":"6","rcsRxGain":"6","rcsLoss":"0","rcsRange":"1.188403","rcsDMin":"0.01","rcsDMax":"20"},"plots":{"fftPlot":{"x0":3369,"x1":3957,"y0":-6.25,"y1":6.127984183175224}},"view":{"start":3369,"stop":3957},"background":{"start":3633,"stop":3692},"analysis":{"start":3704,"stop":3880},"timing":{"period":2.5370000000000003e-05,"framePeriod":0.00317125},"fit":true};
export const initialState=typeof location==='undefined'?null:(decodeState(location.href)??structuredClone(defaultState));
export function restoreControls(root=document) {
  for(const el of root.querySelectorAll('input[id],select[id]')) {
    const v=initialState?.controls[el.id];if(v===undefined)continue;
    if(el.type==='checkbox'){if(typeof v==='boolean')el.checked=v;continue;}
    if(el.tagName==='SELECT'){if([...el.options].some(o=>o.value===String(v)))el.value=String(v);continue;}
    if(typeof v==='string'&&v.length<100 && (el.type!=='number'&&el.type!=='range'||v!==''&&Number.isFinite(Number(v))))el.value=v;
  }
  const db=[...root.querySelectorAll("input[id],select[id]")].find(el=>el.id==="processingLossDb"), legacy=[...root.querySelectorAll("input[id],select[id]")].find(el=>el.id==="phaseLoss");
  if(db && legacy){
    if(initialState?.controls?.processingLossDb===undefined && initialState?.controls?.phaseLoss!==undefined) db.value=String(-10*Math.log10(1-Number(legacy.value)/100));
    legacy.value=String(100*(1-10**(-Number(db.value)/10)));
  }
}
let ready=false,getExtra=()=>({}),timer;
export function saveState() {
  if(!ready)return;
  clearTimeout(timer);timer=setTimeout(()=>{
    const controls=Object.fromEntries([...document.querySelectorAll('input[id],select[id]')].map(el=>[el.id,el.type==='checkbox'?el.checked:el.value]));
    const viewports=Object.fromEntries([...plots].filter(([,p])=>p.viewport).map(([id,p])=>[id,Object.fromEntries(['x0','x1','y0','y1'].map(k=>[k,p.config[k]]))]));
    history.replaceState(null,'',encodeState(location.href,{controls,plots:viewports,...getExtra()}));
  },150);
}
export function initializeState(extra) {
  getExtra=extra;ready=true;
  document.addEventListener('input',saveState);document.addEventListener('change',saveState);
  window.addEventListener('fmcw-view-state',saveState);
  window.addEventListener('fmcw-match',()=>setTimeout(saveState,0));
}
export function registerPlot(plot) {plots.set(plot.root.id,plot);}
export function restorePlot(plot) {
  if(plot.restoredState)return;
  plot.restoredState=true;
  const bounds=initialState?.plots?.[plot.root.id],base=plot.baseConfig;
  if(!bounds||!base)return;
  if(!['x0','x1','y0','y1'].every(k=>Number.isFinite(bounds[k]))||bounds.x0>=bounds.x1||bounds.y0>=bounds.y1)return;
  const clipped={x0:Math.max(base.x0,bounds.x0),x1:Math.min(base.x1,bounds.x1),y0:Math.max(base.y0,bounds.y0),y1:Math.min(base.y1,bounds.y1)};
  if(clipped.x0<clipped.x1&&clipped.y0<clipped.y1)plot.viewport=clipped;
}
