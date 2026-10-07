import * as THREE from 'three';
import {OrbitControls} from './vendor/OrbitControls.js';
import {PROFILES,timing,state,simulate,spectrogram} from './model.mjs';

const $=id=>document.getElementById(id), colours={grid:'#293b4e',text:'#9fb4c9',tx:'#ffca72',q:'#50d0b0',i:'#7facff',total:'#ffe196',beat:'#fa8494',alias:'#c5a3ef',envelope:'#cfdaea'};
let sim,stft,playing=false,fraction=0,last=performance.now(),plots={},pending=false;
const container=$('scene'), scene=new THREE.Scene();scene.background=new THREE.Color('#091525');
const camera=new THREE.PerspectiveCamera(42,1,.01,200), renderer=new THREE.WebGLRenderer({antialias:true});
renderer.setPixelRatio(Math.min(devicePixelRatio,2));container.prepend(renderer.domElement);
const orbit=new OrbitControls(camera,renderer.domElement);orbit.enableDamping=true;orbit.minDistance=.12;orbit.maxDistance=80;
scene.add(new THREE.HemisphereLight(0xb7e0ff,0x182d35,2.4));const light=new THREE.DirectionalLight(0xffffff,2);light.position.set(2,4,2);scene.add(light);
const world=new THREE.Group();scene.add(world);
const radar=new THREE.Mesh(new THREE.BoxGeometry(.10,.12,.07),new THREE.MeshStandardMaterial({color:0x7facff,metalness:.4,roughness:.4}));world.add(radar);
const dish=new THREE.Mesh(new THREE.CylinderGeometry(.04,.03,.015,28),new THREE.MeshStandardMaterial({color:0xb7cdff}));dish.rotation.x=Math.PI/2;radar.add(dish);dish.position.z=-.05;
const target=new THREE.Mesh(new THREE.SphereGeometry(.022,24,16),new THREE.MeshStandardMaterial({color:0xffca72,emissive:0x624413}));world.add(target);
function line(colour,dashed=false){const mat=dashed?new THREE.LineDashedMaterial({color:colour,dashSize:.035,gapSize:.02}):new THREE.LineBasicMaterial({color:colour});const obj=new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(),new THREE.Vector3()]),mat);world.add(obj);return obj;}
const trajectory=line(0xffca72),along=line(0x7facff,true),cross=line(0x7facff,true);
function setLine(obj,a,b){obj.geometry.setFromPoints([a,b]);obj.computeLineDistances();}
let grid,labels=[],fitRadius=1;
function label(text,color) {
  const c=document.createElement('canvas');c.width=512;c.height=80;const x=c.getContext('2d');x.font='34px system-ui';x.fillStyle=color;x.fillText(text,8,48);
  const texture=new THREE.CanvasTexture(c),sprite=new THREE.Sprite(new THREE.SpriteMaterial({map:texture,transparent:true,depthTest:false}));sprite.scale.set(fitRadius*.45,fitRadius*.07,1);world.add(sprite);labels.push(sprite);return sprite;
}
let radarLabel,targetLabel;
function setView() {camera.position.set(orbit.target.x+fitRadius*.8,fitRadius*1.5,orbit.target.z+fitRadius*1.65);orbit.update();}
function geometry() {
  const end=sim.g.v*sim.duration,xmin=Math.min(0,end,sim.g.x0),xmax=Math.max(0,end,sim.g.x0);
  fitRadius=Math.max(.45,Math.hypot(xmax-xmin,sim.g.y0)*.72);
  if(grid){world.remove(grid);grid.geometry.dispose();grid.material.dispose();}
  grid=new THREE.GridHelper(Math.max(2,Math.ceil(fitRadius*3)),20,0x344e65,0x1a2e42);world.add(grid);
  grid.position.set((xmin+xmax)/2,0,sim.g.y0/2);
  radar.position.set(sim.g.x0,.08,sim.g.y0);radar.scale.setScalar(Math.max(1,fitRadius*.35));
  target.scale.setScalar(Math.max(1,fitRadius*.30));
  setLine(trajectory,new THREE.Vector3(-Math.sign(end||1)*Math.max(.07,Math.abs(end)*.06),.025,0),new THREE.Vector3(end+Math.sign(end||1)*Math.max(.07,Math.abs(end)*.06),.025,0));
  setLine(along,new THREE.Vector3(0,.02,0),new THREE.Vector3(sim.g.x0,.02,0));
  setLine(cross,new THREE.Vector3(sim.g.x0,.02,0),new THREE.Vector3(sim.g.x0,.02,sim.g.y0));
  for(const l of labels){world.remove(l);l.material.map.dispose();l.material.dispose();}labels=[];
  radarLabel=label('Radar (x₀, y₀)',colours.i);radarLabel.position.copy(radar.position).add(new THREE.Vector3(0,fitRadius*.12,0));
  targetLabel=label('Target (vt, 0)',colours.tx);
  const xLabel=label('x₀',colours.i);xLabel.position.set(sim.g.x0/2,.06,-fitRadius*.10);
  const yLabel=label('y₀',colours.i);yLabel.position.set(sim.g.x0+fitRadius*.07,.06,sim.g.y0/2);
  const zero=label('0 m',colours.text);zero.position.set(0,.04,-fitRadius*.15);
  orbit.target.set((xmin+xmax)/2,.03,sim.g.y0/2);setView();
}
function setupPlot(id,ylo,yhi,yLabel) {
  const canvas=$(id),ratio=Math.min(devicePixelRatio,2),w=canvas.clientWidth,h=canvas.clientHeight;
  canvas.width=Math.round(w*ratio);canvas.height=Math.round(h*ratio);
  const backing=document.createElement('canvas');backing.width=canvas.width;backing.height=canvas.height;
  const ctx=backing.getContext('2d');ctx.scale(ratio,ratio);
  const L=67,R=15,T=13,B=46,W=w-L-R,H=h-T-B;
  const x=t=>L+t/sim.duration*W,y=v=>T+(yhi-v)/(yhi-ylo)*H;
  ctx.fillStyle='#101d2b';ctx.fillRect(0,0,w,h);ctx.font='10px system-ui';ctx.fillStyle=colours.text;ctx.strokeStyle=colours.grid;ctx.lineWidth=1;
  for(let j=0;j<=4;j++){const val=ylo+(yhi-ylo)*j/4,py=y(val);ctx.beginPath();ctx.moveTo(L,py);ctx.lineTo(w-R,py);ctx.stroke();ctx.textAlign='right';ctx.fillText(Math.abs(val)<.000001?'0':val.toFixed(Math.abs(yhi-ylo)<5?2:1),L-8,py+3);}
  for(let j=0;j<=5;j++){const t=sim.duration*j/5,px=x(t);ctx.beginPath();ctx.moveTo(px,T);ctx.lineTo(px,h-B);ctx.stroke();ctx.textAlign='center';ctx.fillText((t*1e6).toFixed(1),px,h-B+17);}
  ctx.fillText('Elapsed time (µs) · gaps retained',L+W/2,h-7);
  ctx.save();ctx.translate(13,T+H/2);ctx.rotate(-Math.PI/2);ctx.fillText(yLabel,0,0);ctx.restore();
  ctx.save();ctx.beginPath();ctx.rect(L,T,W,H);ctx.clip();
  for(let k=0;k<sim.n;k++){
    ctx.fillStyle='#56677d24';ctx.fillRect(x(k*sim.clock.period+sim.p.adcStart+sim.clock.live),T,(sim.clock.tail??(sim.p.tail+sim.p.idle))/sim.duration*W,H);
    if(sim.n<=5){ctx.fillStyle='#8193a66b';ctx.font='9px system-ui';ctx.textAlign='center';ctx.fillText(`C${k+1}`,x(k*sim.clock.period+sim.clock.live/2),T+11);}
  }
  ctx.restore();
  const p={canvas,backing,ctx,ratio,w,h,L,R,T,B,W,H,x,y,ylo,yhi};plots[id]=p;return p;
}
function curve(p,records,key,color,dashed=false,jump=Infinity) {
  const c=p.ctx;c.save();c.beginPath();c.rect(p.L,p.T,p.W,p.H);c.clip();c.strokeStyle=color;c.lineWidth=1.35;c.setLineDash(dashed?[4,4]:[]);
  for(const record of records){c.beginPath();let previous=null;for(const pt of record){const val=typeof key==='function'?key(pt):pt[key];if(!Number.isFinite(val)){previous=null;continue;}if(previous===null || Math.abs(val-previous)>jump)c.moveTo(p.x(pt.t),p.y(val));else c.lineTo(p.x(pt.t),p.y(val));previous=val;}c.stroke();}c.restore();
}
function drawPlots() {
  const rfTop=(sim.p.f0+sim.p.slope*sim.clock.ramp)/1e9;
  let p=setupPlot('rf',sim.p.f0/1e9-.015,rfTop+.015,'RF frequency (GHz)');
  curve(p,sim.dense,s=>s.tx/1e9,colours.tx);
  curve(p,sim.dense,s=>s.sameRamp?(s.tx+s.frequency)/1e9:NaN,colours.q);
  const maxIF=Math.max(1,...sim.dense.flat().map(s=>Math.abs(s.frequency)/1e6),sim.p.fs/2e6)*1.12;
  p=setupPlot('frequency',-maxIF,maxIF,'Signed frequency (MHz)');
  curve(p,sim.dense,s=>s.frequency/1e6,colours.total);
  curve(p,sim.dense,s=>s.doppler/1e6,colours.i);
  curve(p,sim.dense,s=>s.beat/1e6,colours.beat);
  curve(p,sim.records,s=>s.aliased/1e6,colours.alias,true,sim.p.fs/2e6);
  const maxVoltage=Math.max(...sim.points.map(s=>Math.hypot(s.i,s.q)),1e-6)*1.13;
  p=setupPlot('voltage',-maxVoltage,maxVoltage,'Relative I, Q · return at 1 m = 1');
  curve(p,sim.records,'i',colours.i);curve(p,sim.records,'q',colours.q);
  curve(p,sim.records,s=>Math.hypot(s.i,s.q),colours.envelope,true);
  curve(p,sim.records,s=>-Math.hypot(s.i,s.q),colours.envelope,true);
  p=setupPlot('spectrum',-sim.p.fs/2e6,sim.p.fs/2e6,'Sampled frequency (MHz)');
  const c=p.ctx;c.save();c.beginPath();c.rect(p.L,p.T,p.W,p.H);c.clip();
  for(const frame of stft.frames) {
    const width=stft.hop/sim.p.fs;
    for(let j=0;j<stft.nfft;j++) {
      const db=10*Math.log10(Math.max(frame.power[j],1e-300)/Math.max(stft.peak,1e-300));
      c.fillStyle=colourMap(Math.max(0,Math.min(1,(db+60)/60)));
      const f=(j-stft.nfft/2)*sim.p.fs/stft.nfft/1e6;
      c.fillRect(p.x(frame.t-width/2),p.y(f+sim.p.fs/stft.nfft/1e6),Math.max(1,width/sim.duration*p.W),Math.ceil(p.H/stft.nfft)+.3);
    }
  }
  c.restore();curve(p,sim.records,s=>s.aliased/1e6,colours.total,true,sim.p.fs/2e6);
  $('stftinfo').textContent=`64 samples · 5.12 µs Hann · Fourier scale ${(stft.resolution/1e3).toFixed(1)} kHz`;
  cursor(fraction*sim.duration);
}
function colourMap(t) {
  const stops=[[12,28,46],[32,54,97],[93,69,153],[76,165,167],[156,208,154],[255,226,139]],p=t*(stops.length-1),a=Math.min(stops.length-2,Math.floor(p)),f=p-a;
  return `rgb(${stops[a].map((v,j)=>Math.round(v+(stops[a+1][j]-v)*f)).join(',')})`;
}
function cursor(t) {
  for(const p of Object.values(plots)) {
    const c=p.canvas.getContext('2d');c.setTransform(1,0,0,1,0,0);c.drawImage(p.backing,0,0);c.setTransform(p.ratio,0,0,p.ratio,0,0);
    c.strokeStyle='#edf4f987';c.lineWidth=1;c.beginPath();c.moveTo(p.x(t),p.T);c.lineTo(p.x(t),p.h-p.B);c.stroke();
  }
}
function update() {
  const p=PROFILES[$('profile').value],g={x0:+$('x0').value,y0:+$('y0').value,v:+$('vnumber').value},n=+$('chirps').value;
  $('x0out').textContent=`${g.x0.toFixed(2)} m`;$('y0out').textContent=`${g.y0.toFixed(2)} m`;$('vout').textContent=`${(g.v/1000).toFixed(3)} km/s`;$('chirpsout').textContent=n;
  sim=simulate(p,g,n,$('mode').value);stft=spectrogram(sim);
  $('chirptrack').replaceChildren(...Array.from({length:n},(_,j)=>{const block=document.createElement('div');block.className='chirp-slot';block.textContent=`C${j+1}`;block.title=`Chirp ${j+1}: ADC ${j*sim.clock.period*1e6}–${(j*sim.clock.period+sim.clock.live)*1e6} µs`;return block;}));
  $('timing').innerHTML=`<strong>One chirp</strong><br>ADC: 10.24 µs · 128 complex samples<br>Ramp tail: 0.37 µs · idle: 3.10 µs<br>Repetition: ${(sim.clock.period*1e6).toFixed(2)} µs<br><strong>This sequence</strong><br>ADC live: ${(n*sim.clock.live*1e6).toFixed(2)} µs<br>Through last ADC window: ${(sim.span*1e6).toFixed(2)} µs<br>With final padding: ${(sim.duration*1e6).toFixed(2)} µs`;
  const outside=sim.points.filter(s=>Math.abs(s.frequency)>p.fs/2).length;
  $('status').textContent=sim.mode==='receiver'?'Receiver approximation active · HPFs + range-side mask':outside?'Analogue IF exceeds ±6.25 MHz in some samples; ADC aliasing is shown.':'Ideal mixer · sampled complex I/Q';
  geometry();drawPlots();frameTime();
}
function frameTime() {
  if(!sim)return;const t=fraction*sim.duration,done=fraction>=1,u=done?sim.clock.period-1e-14:t%sim.clock.period,active=!done&&u>=sim.p.adcStart&&u<sim.p.adcStart+sim.clock.live,k=Math.min(sim.n-1,Math.floor(t/sim.clock.period));
  const s=state(t,Math.min(u,sim.clock.ramp),sim.p,sim.g);
  target.position.set(sim.g.v*t,.06,0);targetLabel.position.copy(target.position).add(new THREE.Vector3(0,fitRadius*.10,-fitRadius*.04));
  $('timeout').textContent=`${(t*1e6).toFixed(2)} µs`;
  const liveText=done?'Sequence complete':active?`Chirp ${k+1} / ${sim.n} · ADC LIVE`:u<sim.clock.ramp?`Chirp ${k+1} / ${sim.n} · RAMP TAIL · ADC off`:`IDLE · ADC off${k+1<sim.n?` · next: chirp ${k+2}`:' · end of sequence'}`;
  $('live').textContent=liveText;$('live').classList.toggle('adc-live',active);
  [...$('chirptrack').children].forEach((el,j)=>{el.classList.toggle('current',j===k&&!done);el.classList.toggle('sampling',j===k&&active);el.style.setProperty('--progress',`${j<k?100:j===k?Math.min(100,u/sim.clock.period*100):0}%`);});
  $('geometryread').textContent=`Range now: ${Math.hypot(sim.g.v*t-sim.g.x0,sim.g.y0).toFixed(3)} m · radial velocity: ${(s.vr/1000).toFixed(3)} km/s · chirp ${k+1}/${sim.n} · ${active?'ADC sampling':u<sim.clock.ramp?'ramp tail · ADC off':'idle · ADC off'}`;
  $('time').value=Math.round(fraction*1000);cursor(t);
}
function queue(){if(!pending){pending=true;requestAnimationFrame(()=>{pending=false;update();});}}
for(const id of ['x0','y0','v']) {
  const range=$(id),number=$(id+'number');
  range.addEventListener('input',()=>{number.value=range.value;queue();});
  number.addEventListener('input',()=>{const val=+number.value;if(Number.isFinite(val)&&val>=+range.min&&val<=+range.max){range.value=val;queue();}});
}
$('profile').addEventListener('change',()=>{$('v').value=PROFILES[$('profile').value].speed;$('vnumber').value=PROFILES[$('profile').value].speed;queue();});
for(const id of ['mode','chirps'])$(id).addEventListener('input',queue);
$('play').addEventListener('click',()=>{playing=!playing;$('play').textContent=playing?'Pause':'Play';last=performance.now();});
$('time').addEventListener('input',()=>{fraction=+$('time').value/1000;playing=false;$('play').textContent='Play';frameTime();});
$('reset').addEventListener('click',()=>{for(const [id,value] of Object.entries({x0:.12,y0:.25,v:PROFILES[$('profile').value].speed})){ $(id).value=value;$(id+'number').value=value;}fraction=0;queue();});
$('viewreset').addEventListener('click',setView);
new ResizeObserver(()=>{const w=container.clientWidth,h=container.clientHeight;renderer.setSize(w,h);camera.aspect=w/h;camera.updateProjectionMatrix();if(sim)drawPlots();}).observe(container);
window.addEventListener('resize',()=>{if(sim)drawPlots();});
update();
function animate(now){requestAnimationFrame(animate);if(playing){fraction=(fraction+Math.min(now-last,100)/7000)%1;frameTime();}last=now;orbit.update();renderer.render(scene,camera);}
requestAnimationFrame(animate);
