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
export const initialState=typeof location==='undefined'?null:decodeState(location.href);
export function restoreControls(root=document) {
  for(const el of root.querySelectorAll('input[id],select[id]')) {
    const v=initialState?.controls[el.id];if(v===undefined)continue;
    if(el.type==='checkbox'){if(typeof v==='boolean')el.checked=v;continue;}
    if(el.tagName==='SELECT'){if([...el.options].some(o=>o.value===String(v)))el.value=String(v);continue;}
    if(typeof v==='string'&&v.length<100 && (el.type!=='number'&&el.type!=='range'||v!==''&&Number.isFinite(Number(v))))el.value=v;
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
