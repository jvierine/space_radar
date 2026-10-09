import {restoreControls} from './gui-state.mjs?v=20261009export8';
import {Heatmap} from './plots.mjs?v=20261009export8';
import {estimateRcs,diameterRoots} from './rcs.mjs?v=20261009rcs1';
const section=document.createElement('section');
section.className='card';
section.innerHTML=`<h2>Four-antenna coherent beamforming</h2>
<p>Equal-amplitude phase search + Nelder–Mead; RX0 phase fixed at 0°.</p>
<div class="fields"><label>Phase steps per receiver<input id="beamSteps" type="number" min="2" max="32" step="1" value="10"></label><button id="beamSearch" disabled>Search 1,000 phase combinations</button></div>
<p id="beamStatus" role="status">Run a trajectory fit to compute beamforming and RCS.</p>
<div id="beamResults" hidden><div id="beamTable"></div>
<p id="beamTotal"></p>
<details><summary>Individual receiver signal and quiet matched powers</summary><div id="beamDiagnostics" class="beam-growth"></div><p class="note">Signal: independent RX mean removal. Noise: raw quiet I/Q.</p></details>
<h3>Gain from adding receivers</h3><div id="beamGrowth" class="beam-growth"></div>
<p class="note">Added gain uses separately optimized subsets. Ideal increment = added RX’s linear SNR. Negative = reduced gain; &gt;100% can result from correlated noise.</p>
<div class="beam-marginals">
<div><h3>RX1 × RX2 · MAX over RX3</h3><div id="beamPlot12" class="heatmap compact"></div></div>
<div><h3>RX1 × RX3 · MAX over RX2</h3><div id="beamPlot13" class="heatmap compact"></div></div>
<div><h3>RX2 × RX3 · MAX over RX1</h3><div id="beamPlot23" class="heatmap compact"></div></div>
</div></div>
`;
document.getElementById('matchSettings').after(section);
const rcsSection=document.createElement('div');
rcsSection.innerHTML=`<h3>RCS and equivalent metallic-sphere diameter</h3>
<div class="fields">
<label>System noise (K)<input id="rcsTemperature" type="number" value="9000" min="1"></label>
<label>TX power (dBm)<input id="rcsPower" type="number" value="12" step="0.1"></label>
<label>TX gain (dBi)<input id="rcsTxGain" type="number" value="6" step="0.1"></label>
<label>RX gain per antenna (dBi)<input id="rcsRxGain" type="number" value="6" step="0.1"></label>
<label>Extra loss (dB)<input id="rcsLoss" type="number" value="0" step="0.1"></label>
<label>Distance (m)<input id="rcsRange" type="number" min="0.0001" step="0.001"></label>
<label>Diameter minimum (mm)<input id="rcsDMin" type="number" min="0.001" value="0.01" step="0.01"></label>
<label>Diameter maximum (mm)<input id="rcsDMax" type="number" max="100" value="20" step="0.1"></label>
</div><p id="rcsSummary" class="note"></p><div id="rcsTable" class="beam-growth"></div>
<p class="note">Thermal-noise estimate with ideal 4× RX gain. PEC Mie inversion; multiple diameters are alternative solutions.</p>`;
section.querySelector('#beamResults').insertBefore(rcsSection,section.querySelector('.beam-marginals'));
restoreControls(section);
const button=section.querySelector('#beamSearch'),steps=section.querySelector('#beamSteps'),status=section.querySelector('#beamStatus'),results=section.querySelector('#beamResults');
const pairs=[[1,2,3],[1,3,2],[2,3,1]];
const maps=pairs.map(([x,y])=>new Heatmap(`beamPlot${x}${y}`));
const worker=new Worker('beam-worker.mjs?v=20261009export8',{type:'module'});
let current=null,mainBusy=false,beamBusy=false,id=0,lastBeam=null,pendingAutomatic=false;
function renderRcs() {
 if(!lastBeam || !current)return;
 const value=id=>Number(section.querySelector('#'+id).value),p=current.meta.parameters;
 const frequency=p.f_start+p.freq_slope*(p.T_adc+(current.meta.samples-1)/(2*p.fs));
 const assumptions={temperature:value('rcsTemperature'),range:value('rcsRange'),frequency,time:lastBeam.effectiveTime,txPowerDbm:value('rcsPower'),txGainDbi:value('rcsTxGain'),rxGainDbi:value('rcsRxGain'),lossDb:value('rcsLoss')};
 const minimum=value('rcsDMin')/1000,maximum=value('rcsDMax')/1000;
 const summary=section.querySelector('#rcsSummary'),table=section.querySelector('#rcsTable');
 if(!Object.values(assumptions).every(Number.isFinite) || assumptions.temperature<=0 || assumptions.range<=0 || !(minimum>=1e-6 && maximum>=minimum && maximum<=.1)) {summary.textContent='Enter positive temperature/range and diameter bounds from 0.001 to 100 mm.';table.innerHTML='';return;}
 summary.textContent=`${(frequency/1e9).toFixed(3)} GHz · sampled integration ${(assumptions.time*1e6).toFixed(2)} µs · distance ${assumptions.range.toFixed(4)} m · assumed TX/gains.`;
 const scores=[...lastBeam.single,lastBeam.peak];
 table.innerHTML=`<table><thead><tr><th>Receiver / combination</th><th>Estimated RCS (m²)</th><th>RCS (dBsm)</th><th>Equivalent PEC diameter(s) (mm)</th></tr></thead><tbody>${scores.map((score,i)=>{
  const sigma=estimateRcs(score,assumptions,i===4?4:1),roots=diameterRoots(sigma,frequency,minimum,maximum);
  return `<tr><td>${i<4?'RX'+i:'Four-RX beamformed'}</td><td>${sigma.toExponential(3)}</td><td>${sigma>0?(10*Math.log10(sigma)).toFixed(2):'—'}</td><td>${roots.length?roots.map(d=>(d*1000).toFixed(3)).join(', '):sigma>0?'No root within bounds':'Below noise floor'}</td></tr>`;
 }).join('')}</tbody></table>`;
}
for(const input of rcsSection.querySelectorAll('input'))input.onchange=renderRcs;
const refresh=()=>{button.disabled=!current||mainBusy||beamBusy;steps.disabled=mainBusy||beamBusy;button.textContent=`Search ${(Number(steps.value)**3).toLocaleString()} phase combinations`;};
steps.oninput=refresh;
window.addEventListener('fmcw-match',e=>{id++;current=e.detail;beamBusy=false;pendingAutomatic=true;results.hidden=true;status.textContent=`Ready for ${current.match.pulses} chirps, ${current.match.start}–${current.match.start+current.match.pulses-1}.`;refresh();if(current.match.beam){pendingAutomatic=false;worker.onmessage({data:{id,type:"result",result:current.match.beam}});}else if(!mainBusy)startBeam();});
window.addEventListener('fmcw-invalidated',()=>{id++;current=null;beamBusy=false;pendingAutomatic=false;results.hidden=true;status.textContent='Run a trajectory fit to compute beamforming and RCS.';refresh();});
window.addEventListener('fmcw-busy',e=>{mainBusy=e.detail;refresh();if(!mainBusy && pendingAutomatic)startBeam();});
function startBeam(){
 if(!current || mainBusy || beamBusy)return;
 pendingAutomatic=false;
 const n=Number(steps.value);
 if(!Number.isInteger(n)||n<2||n>32){status.textContent='Choose an integer from 2 to 32 phase steps.';return;}
 beamBusy=true;refresh();worker.postMessage({id:++id,...current,steps:n});
}
button.onclick=startBeam;
worker.onmessage=({data:m})=>{
 if(m.id!==id)return;
 if(m.type==='progress'){status.textContent=m.text;return;}
 beamBusy=false;refresh();
 if(m.type==='error'){status.textContent=m.text;return;}
 const r=m.result,db=v=>10*Math.log10(v),peak=db(r.peak),single=r.single.map(db),best=Math.max(...single);
 lastBeam=r;
 const fit=current.match;
 section.querySelector('#rcsRange').value=(fit.model==='radial-quadratic'?fit.best[0]:Math.hypot(fit.best[2]*fit.midpoint-fit.best[0],fit.best[1])).toFixed(6);
 status.textContent=`RX phases: ${r.phases.map(p=>(p*180/Math.PI).toFixed(1)+'°').join(', ')} · grid ${db(r.gridPeak).toFixed(2)} → refined ${peak.toFixed(2)} dB${r.refinement.converged?'':' · refinement limit'}`;
 section.querySelector('#beamTable').innerHTML=`<table><thead><tr><th>Receiver / combination</th><th>Background-referenced matched energy (dB)</th></tr></thead><tbody>${single.map((v,i)=>`<tr><td>RX${i}</td><td>${v.toFixed(2)}</td></tr>`).join('')}<tr><td>Four-RX phase grid peak</td><td>${db(r.gridPeak).toFixed(2)}</td></tr><tr><td>Nelder–Mead refined four-RX phases</td><td>${peak.toFixed(2)}</td></tr></tbody></table>`;
 const format=(v,d=2)=>Number.isFinite(v)?v.toFixed(d):'—';
 section.querySelector('#beamDiagnostics').innerHTML=`<table><thead><tr><th>Receiver</th><th>Observed matched power (ADC projection²)</th><th>Raw quiet matched power (ADC projection²)</th><th>Observed / RX0 (dB)</th><th>Quiet / RX0 (dB)</th></tr></thead><tbody>${r.channelPower.map((p,i)=>`<tr><td>RX${i}</td><td>${p.observed.toExponential(3)}</td><td>${p.quiet.toExponential(3)}</td><td>${format(db(p.observed/r.channelPower[0].observed))}</td><td>${format(db(p.quiet/r.channelPower[0].quiet))}</td></tr>`).join('')}</tbody></table>`;
 section.querySelector('#beamTotal').textContent=`Total gain: ${format(db(r.total.gain))} dB · ${format(r.total.percentIdeal,1)}% of ideal 6.02 dB. Baseline: ${format(db(r.total.singleAverage))} dB (four-RX linear average).`;
 section.querySelector('#beamGrowth').innerHTML=`<table><thead><tr><th>Receivers</th><th>Matched energy / quiet (dB)</th><th>Added gain (dB)</th><th>Ideal added gain (dB)</th><th>% of ideal increment</th></tr></thead><tbody>${r.growth.map(g=>`<tr><td>${g.count===1?'RX0':`RX0–RX${g.count-1} (+RX${g.count-1})`}</td><td>${format(db(g.score))}</td><td>${g.count>1?format(db(g.gain)):'—'}</td><td>${g.count>1?format(db(g.idealGain)):'—'}</td><td>${g.count>1?format(g.percentIdeal,1):'—'}</td></tr>`).join('')}</tbody></table>`;
 results.hidden=false;
 renderRcs();
 const step=360/r.steps;
 const lo=Math.min(...r.projections.flatMap(p=>[...p].filter(Number.isFinite)));
 maps.forEach((map,i)=>{
  const [x,y,other]=pairs[i];
  map.resetZoom();
  map.set(r.projections[i],r.steps,r.steps,{x0:-step/2,x1:360-step/2,y0:-step/2,y1:360-step/2,xLabel:`RX${x} phase (°; modulo 360)`,yLabel:`RX${y} phase (°; modulo 360)`,topLabel:`MAX over RX${other}`,lo,hi:db(r.gridPeak),kind:1,colorDecimals:2,colorLabel:"Matched energy / quiet (dB)",marker:{column:r.indices[x-1],row:r.indices[y-1]}});
 });
};
refresh();
