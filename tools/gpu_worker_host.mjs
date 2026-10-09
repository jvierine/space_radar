// Node worker adapter for testing the shipped browser DSP worker without a browser.
import fs from 'node:fs';
import {parentPort,workerData} from 'node:worker_threads';
import {pathToFileURL} from 'node:url';
if(workerData.gpu){
 const {create,globals}=await import(pathToFileURL(process.env.WEBGPU_MODULE??'/tmp/fmcw-webgpu-test/node_modules/webgpu/index.js'));
 Object.assign(globalThis,globals);
 Object.defineProperty(globalThis,'navigator',{value:{gpu:create(['backend=metal'])}});
}
globalThis.postMessage=message=>parentPort.postMessage(message);
globalThis.onmessage=null;
globalThis.fetch=async path=>{
 const name=path.split('?')[0];
 return new Response(fs.readFileSync(new URL('../web/lab/'+name,import.meta.url)),{headers:{'Content-Type':name.endsWith('.wasm')?'application/wasm':'application/octet-stream'}});
};
await import('../web/lab/worker.mjs');
parentPort.on('message',data=>globalThis.onmessage({data}));
