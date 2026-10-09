import fs from 'node:fs';
import vm from 'node:vm';
import assert from 'node:assert/strict';
import {estimateRcs,diameterRoots,sphereRcs,snrDb} from '../web/lab/rcs.mjs';
import {LinePlot} from '../web/lab/plots.mjs';
const source=fs.readFileSync('web/lab/app.mjs','utf8');
const body=source.slice(source.indexOf('function updateScan('),source.indexOf("document.addEventListener('change',event=>"));
const values={scanStart:'1000',scanStop:'1010',rcsTemperature:'9000',rcsPower:'12',rcsTxGain:'6',rcsRxGain:'6',rcsLoss:'0',rcsDMin:'.01',rcsDMax:'20'};
const ids=['scanRange','scanVelocity','scanAcceleration','scanPhases','scanSnr','scanRcs','scanDiameter'];
const scanPlots=Object.fromEntries(ids.map(id=>[id,{set(series,config){this.series=series;this.config=config;}}]));
const meta={samples:225,parameters:{f_start:77e9,freq_slope:99e12,T_adc:2e-6,fs:12.5e6}};
const context={$:()=>({disabled:false}),scanPoints:new Map(),scanPlots,meta,document:{getElementById:id=>({value:values[id]})},num:id=>Number(values[id]),clock:k=>(k*25.37e-6).toFixed(6),estimateRcs,diameterRoots,snrDb,db:x=>10*Math.log10(x)};
vm.createContext(context);vm.runInContext(body,context);
const frequency=meta.parameters.f_start+meta.parameters.freq_slope*(meta.parameters.T_adc+224/(2*meta.parameters.fs));
const sigma=sphereRcs(frequency,.004);
for(const start of [1000,1001]){
  const range=2+(start-1000)*.01,time=144e-6;
  const settings={temperature:9000,range,frequency,time,txPowerDbm:12,txGainDbi:6,rxGainDbi:6,lossDb:0};
  const ratio=count=>1+sigma/estimateRcs(2,settings,count);
  context.input={start,midpoint:95e-6,model:'radial-quadratic',best:[range,-300,1e5],beam:{phases:[0,.1,.2,.3],single:Array(4).fill(ratio(1)),peak:ratio(4),effectiveTime:time}};
  vm.runInContext('updateScan(input)',context);
}
for(const id of ['scanRange','scanVelocity','scanAcceleration','scanPhases','scanRcs','scanSnr','scanDiameter'])assert.equal(scanPlots[id].config.mode,'scatter');

for(const series of scanPlots.scanRcs.series){assert.equal(series.y.length,2);for(const v of series.y)assert(Math.abs(v-10*Math.log10(sigma))<1e-8);}
assert(scanPlots.scanDiameter.series.some(series=>series.y.some(v=>Math.abs(v-4)<1e-5)));
assert.equal(scanPlots.scanDiameter.series.length,5*diameterRoots(sigma,frequency,.00001,.02).length,'Every Mie solution is plotted');
values.rcsTemperature='18000';vm.runInContext('drawScanHistory()',context);
assert(Math.abs(scanPlots.scanRcs.series[0].y[0]-10*Math.log10(sigma)-10*Math.log10(2))<1e-8);
assert.equal(context.scanPoints.size,2);
const html=fs.readFileSync('web/lab/index.html','utf8');
assert(html.slice(html.indexOf('id="matchSettings"'),html.indexOf('id="scanTimeline"')).includes('</section>'));
assert(fs.readFileSync('web/lab/beam-ui.mjs','utf8').includes("getElementById('matchSettings').after(section)"));
// Exercise the renderer: scatter series must never connect their points.
let inSeries=false,lines=0,arcs=0;
const ctx=new Proxy({measureText:()=>({width:20}),clip(){inSeries=true;},restore(){inSeries=false;},lineTo(){if(inSeries)lines++;},arc(){arcs++;}}, {get:(object,key)=>key in object?object[key]:()=>{}});
globalThis.devicePixelRatio=1;
const plot=Object.create(LinePlot.prototype);
plot.root={clientWidth:1200,clientHeight:220};plot.canvas={getContext:()=>ctx};
plot.config={x0:0,x1:2,mode:"scatter"};plot.series=[{name:'test',color:'#000',x:[0,1,2],y:[1,2,NaN]}];
plot.draw();assert.equal(lines,0);assert.equal(arcs,2);
plot.config.mode='line';lines=0;arcs=0;plot.draw();assert.equal(lines,1);assert.equal(arcs,0);
console.log('PASS: unconnected scatter rendering, five RCS/diameter histories, every Mie root, beam gain, recalibration and section ordering');
