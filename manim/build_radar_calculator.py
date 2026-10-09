"""Generate a self-contained parameter calculator from verified HDF5 numbers."""
from pathlib import Path
import json
import h5py
HERE=Path(__file__).resolve().parent
with h5py.File(HERE/'assets/fmcw_lecture.h5') as h:
    models=[dict(name='3 mm sphere, 78 GHz',snr=h['three_mm/snr_gain_6_db'][0,0].item()),dict(name='6606: 5 mm, 78 GHz',snr=h['parallel/ideal_1m_snr_db'][0].item()),dict(name='6607: 6 mm, 79.5 GHz',snr=h['perpendicular/ideal_1m_snr_db'][0].item())]
page='''<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width"><meta name="author" content="Juha Vierinen"><title>Space radar SNR calculator</title>
<style>body{background:#07111f;color:#f3f7fa;font:20px Georgia;margin:3em auto;max-width:1000px;padding:0 1em}h1{color:#f3d35a}label{display:block;margin:1.2em 0}input{width:280px;margin:0 1em}select{font:inherit}table{border-collapse:collapse;width:100%;margin:1.5em 0}th,td{padding:.8em;text-align:right;border-bottom:1px solid #456}th:first-child,td:first-child{text-align:left}small{color:#9db0c3}a{color:#74a9ff}output{color:#45c2b1}</style>
<h1>Space radar: ideal thermal SNR</h1><p>Juha Vierinen</p><p>Change transmitted power and antenna gain. The table uses coherent integration of 1, 2, 4 or 10 recorded chirps.</p>
<label>Target <select id="target"></select></label>
<label>TX power <input id="pt" type="range" min="-10" max="12" step=".5" value="12"><output id="ptv"></output></label>
<label>TX gain <input id="gt" type="range" min="0" max="20" step=".5" value="6"><output id="gtv"></output></label>
<label>RX gain <input id="gr" type="range" min="0" max="20" step=".5" value="6"><output id="grv"></output></label>
<label>System temperature <input id="temp" type="range" min="1000" max="20000" step="100" value="10000"><output id="tempv"></output></label>
<label>Extra loss <input id="loss" type="range" min="0" max="20" step=".5" value="0"><output id="lossv"></output></label>
<table><thead><tr><th>Chirps</th><th>ADC live time</th><th>1 m</th><th>10 m</th><th>100 m</th></tr></thead><tbody id="table"></tbody></table>
<p>At 6.70 km/s a one-metre window lasts 149.3 µs, about 11 chirp periods; at 6.285 km/s, 159.1 µs, about 12 periods. Entry and exit chirps can be partial. Both recorded shots fit a ten-chirp window: approximately 133.63 µs elapsed span and 102.4 µs of ADC data. At 14 km/s the same one-metre window gives only about 5 chirp periods.</p>
<small>One TX, one RX; perfectly conducting sphere; fixed range; perfect waveform match; 128 samples at 12.5 MS/s per chirp. These are conditional thermal budgets, not measured calibrated SNR. Changing chamber echoes, pointing, filtering and range migration affect practical detection. Gain above the board peak would require a different antenna. No extra gain from unobserved gaps. Source: build_radar_calculator.py and space_radar_assets.py / assets/fmcw_lecture.h5.</small>
<p><a href="index.html">Animated lecture</a></p>
<script>const models=MODELS;const get=id=>document.getElementById(id);models.forEach((m,i)=>get('target').add(new Option(m.name,i)));function update(){['pt','gt','gr','temp','loss'].forEach(k=>get(k+'v').textContent=get(k).value+({pt:' dBm',gt:' dBi',gr:' dBi',temp:' K',loss:' dB'}[k]));let base=models[+get('target').value].snr+(+get('pt').value-12)+(+get('gt').value-6)+(+get('gr').value-6)-10*Math.log10(+get('temp').value/10000)-(+get('loss').value);get('table').innerHTML=[1,2,4,10].map(n=>'<tr><td>'+n+'</td><td>'+(n*10.24).toFixed(2)+' µs</td>'+[1,10,100].map(r=>'<td>'+(base+10*Math.log10(n)-40*Math.log10(r)).toFixed(1)+' dB</td>').join('')+'</tr>').join('')}document.querySelectorAll('input,select').forEach(x=>x.addEventListener('input',update));update();</script>'''.replace('MODELS',json.dumps(models))
out=HERE/'qa_fmcw_web/fmcw-space-radar';out.mkdir(parents=True,exist_ok=True);(out/'calculator.html').write_text(page)
print(out/'calculator.html')
index=out/'index.html'
if index.exists():
    text=index.read_text()
    if 'id="radar-calculator-link"' not in text:
        link='<a id="radar-calculator-link" href="calculator.html" target="_blank" style="position:fixed;bottom:8px;left:12px;z-index:9999;color:#9db0c3;font:13px Georgia;text-decoration:none">SNR calculator</a>'
        text=text.replace('</body>',link+'\n</body>')
    if 'id="radar-simulator-link"' not in text:
        link='<a id="radar-simulator-link" href="https://juha.no/fmcw/sim/" target="_blank" rel="noopener" style="position:fixed;bottom:8px;left:135px;z-index:9999;color:#9db0c3;font:13px Georgia;text-decoration:none">FMCW simulator ↗</a>'
        text=text.replace('</body>',link+'\n</body>')
    index.write_text(text)
