// HDF5 schema: every completed train estimate and beam result; current full MF bank.
export async function writeAnalysis(h5,payload) {
  const {FS}=await h5.ready, filename='/fmcw-analysis.h5';
  const file=new h5.File(filename,'w');
  const attr=(g,k,v)=>g.create_attribute(k,typeof v==='string'?v:JSON.stringify(v));
  const dataset=(g,name,data,shape,dtype='<d')=>{
    if(!data?.length)return;
    return g.create_dataset({name,data: dtype==='<f'?Float32Array.from(data):Float64Array.from(data),shape:shape??[data.length],dtype});
  };
  try {
    attr(file,'schema','fmcw-analysis-1');attr(file,'created_utc',new Date().toISOString());
    attr(file,'source_url',payload.url);attr(file,'recording',payload.meta);attr(file,'settings',payload.settings);
    attr(file,'timing','Reconstructed seconds since file start; not measured recording timing.');
    attr(file,'noise_reference','Raw full-sample-bandwidth E[|I+iQ|^2] and receiver covariance; matched noise power = sample noise power times template energy. Events independently mean-subtracted per RX.');
    attr(file,'coverage','All completed train estimates and beamforming results. Full matched-filter bank and voltage fit for current displayed train.');
    const history=file.create_group('history'),points=payload.points,N=points.length;
    if(N){
      dataset(history,'start_chirp',points.map(p=>p.start));dataset(history,'pulses',points.map(p=>p.pulses));
      dataset(history,'time_since_file_start_s',points.map(p=>p.time));
      dataset(history,'r0_v0_a0',points.flatMap(p=>p.kinematics),[N,3]);
      attr(history,'r0_v0_a0_units',['m','m/s','m/s^2']);
      dataset(history,'single_rx_matched_to_quiet',points.flatMap(p=>Array.from(p.beam.single)),[N,4]);
      dataset(history,'beam_matched_to_quiet',points.map(p=>p.beam.peak));
      dataset(history,'phases_rad',points.flatMap(p=>Array.from(p.beam.phases)),[N,4]);
      dataset(history,'rcs_m2',points.flatMap(p=>p.rcs.map(r=>r.sigma)),[N,5]);
      attr(history,'rcs_channel_order',['RX0','RX1','RX2','RX3','four_receiver_beam']);
      const diameters=[],index=[];
      for(const p of points)for(const r of p.rcs){index.push(diameters.length,r.diameters.length);diameters.push(...r.diameters);}
      dataset(history,'diameter_mm',diameters);dataset(history,'diameter_start_count',index,[N,5,2]);
      attr(history,'diameter_model','All PEC metallic-sphere Mie inversion roots within configured bounds.');
      const trains=file.create_group('trains');
      for(const p of points){
        const g=trains.create_group(`chirp_${p.start}`);attr(g,'best',p.best);attr(g,'model',p.model??'radial');attr(g,'grid',p.grid);
        attr(g,'reference_starts',p.referenceStarts??[]);
        // Preserve every beam output, including the three MAX projections, powers and gain increments.
        for(const [key,value] of Object.entries(p.beam)){
          if(Array.isArray(value)&&value.length&&value.every(x=>ArrayBuffer.isView(x)))dataset(g,key,value.flatMap(x=>Array.from(x)),[value.length,value[0].length]);
          else if(ArrayBuffer.isView(value)||Array.isArray(value)&&value.every(x=>typeof x==='number'))dataset(g,key,value);
          else attr(g,key,value);
        }
      }
    }
    if(payload.result){
      const r=payload.result,g=file.create_group('current_match');
      for(const key of ['start','pulses','period','midpoint','model','grid','radialSpec','backend','seconds','searchReceivers','noiseModel','noiseSamples','noisePower'])if(r[key]!==undefined)attr(g,key,r[key]);
      for(const key of ['best','cube','nodes','observed','fit','scanResults'])if(r[key])dataset(g,key,r[key],key==='observed'||key==='fit'?[r[key].length/2,2]:undefined,key==='cube'?'<f':'<d');
      if(r.fit&&r.observed)dataset(g,'residual',Array.from(r.observed,(v,i)=>v-r.fit[i]),[r.fit.length/2,2]);
      attr(g,'cube_order','Flattened (range, acceleration, velocity) for radial search; use radialSpec dimensions and grid bounds.');
    }
  }finally{file.close();}
  const bytes=FS.readFile(filename);FS.unlink(filename);return bytes;
}
