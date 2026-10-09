// MathJax SVG labels are drawn into the axes canvas, including PNG exports.
const cache=new Map();
export function ensureMath(tex) {
  if(!tex)return Promise.resolve(null);
  if(cache.has(tex))return cache.get(tex).promise;
  const entry={image:null,width:0,height:0};cache.set(tex,entry);
  entry.promise=(async()=>{
    if(!globalThis.MathJax?.startup?.promise)return null;
    await MathJax.startup.promise;
    const node=await MathJax.tex2svgPromise(tex,{display:false}),svg=node.querySelector('svg');
    const box=svg.getAttribute('viewBox').split(/\s+/).map(Number);
    const height=parseFloat(svg.getAttribute('height'))*7;
    entry.height=height;entry.width=height*box[2]/box[3];
    svg.setAttribute('width',String(entry.width));svg.setAttribute('height',String(height));
    svg.setAttribute('color','#20252b');svg.style.color='#20252b';
    const image=new Image();
    await new Promise((resolve,reject)=>{image.onload=resolve;image.onerror=reject;image.src='data:image/svg+xml;charset=utf-8,'+encodeURIComponent(new XMLSerializer().serializeToString(svg));});
    entry.image=image;return entry;
  })().catch(()=>null);
  return entry.promise;
}
export function drawMath(plot,ctx,tex,x,y,rotate=false) {
  if(!tex)return false;
  const entry=cache.get(tex);
  if(!entry){ensureMath(tex).then(()=>plot.draw());return false;}
  if(!entry.image)return false;
  ctx.save();ctx.translate(x,y);if(rotate)ctx.rotate(-Math.PI/2);
  ctx.drawImage(entry.image,-entry.width/2,-entry.height/2,entry.width,entry.height);ctx.restore();return true;
}
