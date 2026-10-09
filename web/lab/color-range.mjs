// Match each colorbar to finite texels in the displayed viewport.
export function visibleColorRange(values,width,height,base,view) {
  const edge=(v,lo,hi,n,round)=>Math.max(0,Math.min(n,round((v-lo)/(hi-lo)*n)));
  const x0=edge(view.x0,base.x0,base.x1,width,Math.floor);
  const x1=edge(view.x1,base.x0,base.x1,width,Math.ceil);
  const y0=edge(view.y0,base.y0,base.y1,height,Math.floor);
  const y1=edge(view.y1,base.y0,base.y1,height,Math.ceil);
  let lo=Infinity,hi=-Infinity;
  for(let y=y0;y<y1;y++)for(let x=x0;x<x1;x++) {
    const v=values[y*width+x];
    if(Number.isFinite(v)){lo=Math.min(lo,v);hi=Math.max(hi,v);}
  }
  if(!Number.isFinite(lo))return {lo:0,hi:1};
  return {lo,hi:hi>lo?hi:lo+1};
}
