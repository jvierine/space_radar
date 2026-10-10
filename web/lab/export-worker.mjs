import * as h5 from './vendor/h5wasm/hdf5_hl.js';
import {writeAnalysis} from './export-analysis.mjs?v=20261010range18';
self.onmessage=async({data})=>{
  try{const bytes=await writeAnalysis(h5,data);self.postMessage({bytes},[bytes.buffer]);}
  catch(error){self.postMessage({error:error.message});}
};
