import {restoreControls} from './gui-state.mjs?v=20261010scale21';
import {estimateRcs,diameterRoots,snrDb} from './rcs.mjs?v=20261010loss2db17';
const section=document.createElement('section');
section.className='card';
section.innerHTML=`<h2>Four-antenna coherent beamforming</h2>
<p>Joint noise-weighted least squares: trajectory and four complex RX amplitudes. The same fit supplies the coherent beam weights.</p>
<div class="fields"><input id="beamSteps" type="hidden" value="10"><button id="beamSearch" disabled>Recalculate fitted beam</button></div>
<p id="beamStatus" role="status">Run a trajectory fit to compute beamforming and RCS.</p>
<div id="beamResults" hidden><div id="beamTable"></div>
<p id="beamTotal"></p>
<details><summary>Individual receiver signal and quiet matched powers</summary><div id="beamDiagnostics" class="beam-growth"></div><p class="note">Signal: independent RX mean removal. Noise: raw quiet I/Q.</p></details>
<h3>Gain from adding receivers</h3><div id="beamGrowth" class="beam-growth"></div>
<p class="note">Each subset uses its background covariance. Fitted-amplitude noise power is subtracted; trajectory search selection can still bias weak peaks.</p>
<div class="beam-marginals" hidden></div></div>
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
<p class="note">Thermal-noise estimate with covariance-calibrated combining gain, assuming equal calibrated antenna responses. PEC Mie inversion; multiple diameters are alternative solutions. Min–max is a candidate envelope within the selected bounds, not a confidence interval or a continuous set of solutions.</p>`;
section.querySelector('#beamResults').insertBefore(rcsSection,section.querySelector('.beam-marginals'));
restoreControls(section);
const button=section.querySelector('#beamSearch'),steps=section.querySelector('#beamSteps'),status=section.querySelector('#beamStatus'),results=section.querySelector('#beamResults');
const worker=new Worker('beam-worker.mjs?v=20261010scale21',{type:'module'});
let current=null,mainBusy=false,beamBusy=false,id=0,lastBeam=null,pendingAutomatic=false;
function renderRcs() {
 if(!lastBeam || !current)return;
 const value=id=>Number(section.querySelector('#'+id).value),p=current.meta.parameters;
 const frequency=p.f_start+p.freq_slope*(p.T_adc+(current.meta.samples-1)/(2*p.fs));
 const assumptions={temperature:value('rcsTemperature'),range:value('rcsRange'),frequency,T_coh:lastBeam.T_coh??lastBeam.effectiveTime,txPowerDbm:value('rcsPower'),txGainDbi:value('rcsTxGain'),rxGainDbi:value('rcsRxGain'),lossDb:value('rcsLoss')};
 const minimum=value('rcsDMin')/1000,maximum=value('rcsDMax')/1000;
 const summary=section.querySelector('#rcsSummary'),table=section.querySelector('#rcsTable');
 if(lastBeam.noiseCalibrated===false){summary.textContent='RCS needs at least two background chirps to separate noise from stationary echoes.';table.innerHTML='';return;}
 if(!Object.values(assumptions).every(Number.isFinite) || assumptions.temperature<=0 || assumptions.range<=0 || !(minimum>=1e-6 && maximum>=minimum && maximum<=.1)) {summary.textContent='Enter positive temperature/range and diameter bounds from 0.001 to 100 mm.';table.innerHTML='';return;}
 summary.textContent=`${(frequency/1e9).toFixed(3)} GHz · T_coh ${(assumptions.T_coh*1e6).toFixed(2)} µs · distance ${assumptions.range.toFixed(4)} m · assumed TX/gains.`;
 const scores=[...lastBeam.single,lastBeam.peak];
 table.innerHTML=`<table><thead><tr><th>Receiver / combination</th><th>Estimated RCS (m²)</th><th>RCS (dBsm)</th><th>Candidate envelope (mm)</th><th>Equivalent PEC diameter(s) (mm)</th></tr></thead><tbody>${scores.map((score,i)=>{
  const sigma=estimateRcs(score,assumptions,i===4?lastBeam.rcsGain:1),roots=diameterRoots(sigma,frequency,minimum,maximum);
  return `<tr><td>${i<4?'RX'+i:'Four-RX beamformed'}</td><td>${sigma.toExponential(3)}</td><td>${sigma>0?(10*Math.log10(sigma)).toFixed(2):'—'}</td><td>${roots.length?`${(1000*Math.min(...roots)).toFixed(3)}–${(1000*Math.max(...roots)).toFixed(3)}`:'—'}</td><td>${roots.length?roots.map(d=>(d*1000).toFixed(3)).join(', '):sigma>0?'No root within bounds':'Below noise floor'}</td></tr>`;
 }).join('')}</tbody></table>`;
}
for(const input of rcsSection.querySelectorAll('input'))input.onchange=renderRcs;
const refresh=()=>{button.disabled=!current||mainBusy||beamBusy;steps.disabled=mainBusy||beamBusy;button.textContent="Recalculate fitted beam";};
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
 const r=m.result,db=v=>10*Math.log10(v),peak=snrDb(r.peak),single=r.single.map(snrDb),best=Math.max(...single);
 lastBeam=r;
 const fit=current.match;
 section.querySelector('#rcsRange').value=(fit.model==='radial-quadratic'?fit.best[0]:Math.hypot(fit.best[2]*fit.midpoint-fit.best[0],fit.best[1])).toFixed(6);
 status.textContent=`T_coh ${(r.T_coh*1e6).toFixed(2)} µs · Analysis bandwidth ${(r.analysisBandwidth/1000).toFixed(2)} kHz · Fitted RX phases relative to RX0: ${r.echoPhases.map(p=>(p*180/Math.PI).toFixed(1)+'°').join(', ')}${r.noiseCalibrated===false?' · single-chirp noise upper bound':''}`;
 section.querySelector('#beamTable').innerHTML=`<table><thead><tr><th>Receiver / combination</th><th>SNR in analysis bandwidth (dB)</th></tr></thead><tbody>${single.map((v,i)=>`<tr><td>RX${i}</td><td>${v.toFixed(2)}</td></tr>`).join('')}<tr><td>Joint noise-weighted four-RX fit</td><td>${peak.toFixed(2)}</td></tr></tbody></table>`;
 const format=(v,d=2)=>Number.isFinite(v)?v.toFixed(d):'—';
 section.querySelector('#beamDiagnostics').innerHTML=`<table><thead><tr><th>Receiver</th><th>Background noise power (ADC²)</th><th>Observed power in analysis bandwidth (ADC²)</th><th>Noise in analysis bandwidth (ADC²)</th><th>Channel SNR (dB)</th></tr></thead><tbody>${r.channelPower.map((p,i)=>`<tr><td>RX${i}</td><td>${r.rawNoisePower?.[i]?.toExponential(3)??'—'}</td><td>${r.analysisObservedPower?.[i]?.toExponential(3)??'—'}</td><td>${r.analysisNoisePower?.[i]?.toExponential(3)??'—'}</td><td>${format(snrDb(r.single[i]))}</td></tr>`).join('')}</tbody></table>`;
 section.querySelector('#beamTotal').textContent=`GLS matched score ${format(r.glsScore)} · fitted-amplitude noise contribution ${format(r.fittedNoiseBias)} · covariance ridge ${r.regularization} · calibration gain ${format(r.rcsGain)}. Complex amplitudes and combining weights are solved analytically at the refined trajectory.`;
 section.querySelector('#beamGrowth').innerHTML=`<table><thead><tr><th>Receivers</th><th>Fitted signal SNR (dB)</th></tr></thead><tbody>${r.growth.map(g=>`<tr><td>RX0–RX${g.count-1}</td><td>${format(snrDb(g.score))}</td></tr>`).join('')}</tbody></table>`;
 results.hidden=false;
 renderRcs();

};
refresh();
