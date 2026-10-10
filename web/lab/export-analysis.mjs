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
    attr(file,'schema','fmcw-analysis-2');attr(file,'joint_fit','Bounded Nelder-Mead variable-projection GLS; complex receive amplitudes solved analytically at each motion trial. White temporal noise; full background receive covariance. Beam excess SNR subtracts fitted-amplitude noise contribution; trajectory-search selection bias remains.');;attr(file,'created_utc',new Date().toISOString());
    attr(file,'source_url',payload.url);attr(file,'recording',payload.meta);attr(file,'settings',payload.settings);
    attr(file,'timing','Reconstructed seconds since file start; not measured recording timing.');
    attr(file,'noise_reference','Full-sample-bandwidth per-RX background residual variance and receiver covariance, with M/(M-1) mean-fit correction; single-chirp fallback is an unseparated upper bound; matched noise power = sample noise power times template energy. Events independently mean-subtracted per RX.');
    attr(file,'coverage','All completed train estimates and beamforming results. Full matched-filter bank and voltage fit for current displayed train.');
    const history=file.create_group('history'),points=payload.points,N=points.length;
    if(N){
      dataset(history,'start_chirp',points.map(p=>p.start));dataset(history,'pulses',points.map(p=>p.pulses));
      dataset(history,'time_since_file_start_s',points.map(p=>p.time));
      dataset(history,'r0_v0_a0',points.flatMap(p=>p.kinematics),[N,3]);
      attr(history,'r0_v0_a0_units',['m','m/s','m/s^2']);
      dataset(history,'single_rx_matched_to_quiet',points.flatMap(p=>Array.from(p.beam.single)),[N,4]);
      dataset(history,'T_coh_s',points.map(p=>p.beam.T_coh??p.beam.effectiveTime));
      dataset(history,'analysis_observed_power_adc2',points.flatMap(p=>Array.from(p.beam.analysisObservedPower??[NaN,NaN,NaN,NaN])),[N,4]);
      dataset(history,'analysis_bandwidth_hz',points.map(p=>p.beam.analysisBandwidth??NaN));
      dataset(history,'analysis_noise_power_adc2',points.flatMap(p=>Array.from(p.beam.analysisNoisePower??[NaN,NaN,NaN,NaN])),[N,4]);
      dataset(history,'single_rx_snr_linear',points.flatMap(p=>Array.from(p.beam.single,v=>Math.max(v-1,0))),[N,4]);
      dataset(history,'beam_snr_linear',points.map(p=>Math.max(p.beam.peak-1,0)));
      dataset(history,'raw_noise_power_adc2',points.flatMap(p=>Array.from(p.beam.rawNoisePower??[NaN,NaN,NaN,NaN])),[N,4]);
      dataset(history,'matched_noise_power',points.flatMap(p=>p.beam.channelPower?.map(c=>c.quiet)??[NaN,NaN,NaN,NaN]),[N,4]);
      dataset(history,'beam_matched_to_quiet',points.map(p=>p.beam.peak));
      dataset(history,'gls_matched_score',points.map(p=>p.beam.glsScore??NaN));
      dataset(history,'fitted_amplitude_noise_bias',points.map(p=>p.beam.fittedNoiseBias??NaN));
      dataset(history,'combining_calibration_gain',points.map(p=>p.beam.rcsGain??NaN));
      dataset(history,'phases_rad',points.flatMap(p=>Array.from(p.beam.phases)),[N,4]);
      dataset(history,'rcs_m2',points.flatMap(p=>p.rcs.map(r=>r.sigma)),[N,5]);
      attr(history,'rcs_channel_order',['RX0','RX1','RX2','RX3','four_receiver_beam']);
      const diameters=[],index=[];
      for(const p of points)for(const r of p.rcs){index.push(diameters.length,r.diameters.length);diameters.push(...r.diameters);}
      dataset(history,'diameter_mm',diameters);dataset(history,'diameter_start_count',index,[N,5,2]);
      dataset(history,'diameter_candidate_envelope_mm',points.flatMap(p=>p.rcs.flatMap(r=>r.diameters.length?[Math.min(...r.diameters),Math.max(...r.diameters)]:[NaN,NaN])),[N,5,2]);
      attr(history,'diameter_envelope','Min and max inversion roots within configured diameter bounds; not a confidence interval. Discrete alternative roots retained.');
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
      for(const key of ['start','pulses','period','midpoint','model','grid','radialSpec','backend','seconds','searchReceivers','noiseModel','noiseSamples','noisePower','noiseSamplesPerRx','noiseCalibrated'])if(r[key]!==undefined)attr(g,key,r[key]);
      for(const key of ['best','cube','nodes','observed','fit','scanResults','aliasCandidates'])if(r[key])dataset(g,key,r[key],key==='observed'||key==='fit'?[r[key].length/2,2]:undefined,key==='cube'?'<f':'<d');
      if(r.fit&&r.observed)dataset(g,'residual',Array.from(r.observed,(v,i)=>v-r.fit[i]),[r.fit.length/2,2]);
      attr(g,'cube_order','Flattened (range, acceleration, velocity) for radial search; use radialSpec dimensions and grid bounds.');
    }
  }finally{file.close();}
  const bytes=FS.readFile(filename);FS.unlink(filename);return bytes;
}
