import {Heatmap} from './plots.mjs?v=20261009k';
const section=document.createElement('section');
section.className='card';
section.innerHTML=`<h2>Four-antenna coherent beamforming</h2>
<p>Fix RX0 at 0°. Search the three relative receiver phases with equal amplitudes, after fitting the trajectory. Rotate stored I/Q by exp(+iφ)/2 before combining.</p>
<div class="fields"><label>Phase steps per receiver<input id="beamSteps" type="number" min="2" max="32" step="1" value="10"></label><button id="beamSearch" disabled>Search 1,000 phase combinations</button></div>
<p id="beamStatus" role="status">Fit a trajectory first, then optimize the receiver phases for that train.</p>
<div id="beamResults" hidden><div id="beamTable"></div><h3>RX1 × RX2 phase · maximum over RX3 phase</h3><div id="beamPlot" class="heatmap"></div></div>
<p class="note">This is a conditional phase search at the fitted trajectory, not a joint trajectory search or an arrival-angle estimate. Each receiver has its own complex quiet mean removed using common intact chirps. All combinations use the same independent quiet-reference trains and retain cross-receiver noise correlation. Maximizing 1,000 trials can raise a noise peak: gain here is not independent detection evidence. The measurement plots and trajectory fit above remain the selected single receiver.</p>`;
document.getElementById('matchSettings').after(section);
const button=section.querySelector('#beamSearch'),steps=section.querySelector('#beamSteps'),status=section.querySelector('#beamStatus'),results=section.querySelector('#beamResults');
const map=new Heatmap('beamPlot');
const worker=new Worker('beam-worker.mjs?v=20261009beam2',{type:'module'});
let current=null,mainBusy=false,beamBusy=false,id=0;
const refresh=()=>{button.disabled=!current||mainBusy||beamBusy;steps.disabled=mainBusy||beamBusy;button.textContent=`Search ${(Number(steps.value)**3).toLocaleString()} phase combinations`;};
steps.oninput=refresh;
window.addEventListener('fmcw-match',e=>{id++;current=e.detail;beamBusy=false;results.hidden=true;status.textContent=`Ready for ${current.match.pulses} chirps, ${current.match.start}–${current.match.start+current.match.pulses-1}.`;refresh();});
window.addEventListener('fmcw-invalidated',()=>{id++;current=null;beamBusy=false;results.hidden=true;status.textContent='Fit a trajectory first, then optimize the receiver phases for that train.';refresh();});
window.addEventListener('fmcw-busy',e=>{mainBusy=e.detail;refresh();});
button.onclick=()=>{
 const n=Number(steps.value);
 if(!Number.isInteger(n)||n<2||n>32){status.textContent='Choose an integer from 2 to 32 phase steps.';return;}
 beamBusy=true;refresh();worker.postMessage({id:++id,...current,steps:n});
};
worker.onmessage=({data:m})=>{
 if(m.id!==id)return;
 if(m.type==='progress'){status.textContent=m.text;return;}
 beamBusy=false;refresh();
 if(m.type==='error'){status.textContent=m.text;return;}
 const r=m.result,db=v=>10*Math.log10(v),peak=db(r.peak),single=r.single.map(db),best=Math.max(...single);
 status.textContent=`Best relative phases RX0–RX3: ${r.phases.map(p=>(p*180/Math.PI).toFixed(1)+'°').join(', ')} · peak ${peak.toFixed(2)} dB · ${(peak-best).toFixed(2)} dB versus best individual receiver · ${r.steps**3} phase combinations · ${r.referenceStarts.length} common noise trains · ${r.meanCount} common mean chirps.`;
 section.querySelector('#beamTable').innerHTML=`<table><thead><tr><th>Receiver / combination</th><th>Background-referenced matched energy (dB)</th></tr></thead><tbody>${single.map((v,i)=>`<tr><td>RX${i}</td><td>${v.toFixed(2)}</td></tr>`).join('')}<tr><td>Best four-RX phase combination</td><td>${peak.toFixed(2)}</td></tr></tbody></table>`;
 results.hidden=false;
 const step=360/r.steps;
 map.set(r.projection,r.steps,r.steps,{x0:-step/2,x1:360-step/2,y0:-step/2,y1:360-step/2,xLabel:'RX1 phase (degrees; modulo 360)',yLabel:'RX2 phase (degrees; modulo 360)',topLabel:`Maximum over RX3 · peak ${peak.toFixed(2)} dB`,lo:Math.min(...r.projection.filter(Number.isFinite)),hi:peak,kind:1,colorDecimals:2,colorLabel:"Matched energy / quiet noise (dB)",marker:{column:r.indices[0],row:r.indices[1]}});
};
refresh();
