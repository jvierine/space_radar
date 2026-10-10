export function selectDataset(catalog,url) {
  if(catalog.version!==1 || !Array.isArray(catalog.datasets) || !catalog.datasets.length)throw Error('Recording catalog is empty.');
  const ids=new Set();
  for(const d of catalog.datasets){
    if(!/^[a-zA-Z0-9_-]+$/.test(d.id)||ids.has(d.id)||typeof d.label!=='string'||!/^datasets\/[a-zA-Z0-9_-]+\/$/.test(d.path))throw Error('Invalid recording catalog.');
    ids.add(d.id);
  }
  const id=new URL(url).searchParams.get('dataset')??catalog.datasets[0].id;
  const entry=catalog.datasets.find(d=>d.id===id);
  if(!entry)throw Error(`Unknown recording: ${id}`);
  return entry;
}
export function datasetUrl(url,id){
  const u=new URL(url);u.searchParams.set('dataset',id);u.searchParams.delete('gui');return u.href;
}
export function recordingDefaults(meta){
  const total=meta.total_chirps, per=meta.chirps_per_frame;
  const start=Math.max(0,Math.min(total-per,Math.floor(meta.frame_of_interest??0)*per));
  return {start:0,stop:total,chirp:start,bgStart:0,bgStop:Math.min(per,total),scanStart:start,scanStop:Math.min(total,start+per)};
}

export function voltageColorRange(values){
  let peak=0;
  for(const value of values)if(Number.isFinite(value))peak=Math.max(peak,Math.abs(value));
  return Math.max(1,Math.ceil(peak*1.02));
}
