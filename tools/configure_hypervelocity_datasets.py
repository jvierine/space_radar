"""Configure packaged four-frame excerpts of Ellingsen's shots 6606/6607.
Run from the repository root after prepare_lab_data.py for frames 43:47/34:38.
Quiet windows follow the previously inspected pre-disturbance intervals;
analysis starts are display selections, not independently verified triggers.
"""
import json, hashlib
from pathlib import Path
for shot,onset,quiet in [(6606,3386,3386),(6607,1624,1624)]:
 p=Path(f'web/lab/datasets/shot{shot}/metadata.json')
 m=json.loads(p.read_text())
 m['gui_defaults_generator']='tools/configure_hypervelocity_datasets.py'
 m['gui_defaults_generator_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 m['gui_defaults']={'controls':{'rx':'0','chirp':str(onset),'chirpSlider':str(onset),'algorithm':'fft','pulses':'8','processingLossDb':'2','vMin':'-7000','vMax':'7000','xMin':'0.001','xMax':'3','yMin':'0','yMax':'1000000'},'view':{'start':max(0,onset-500),'stop':min(m['total_chirps'],onset+800)},'background':{'start':0,'stop':quiet},'analysis':{'start':onset,'stop':onset+128}}
 p.write_text(json.dumps(m,indent=2)+'\n')
p=Path('web/lab/dataset-catalog.json');c=json.loads(p.read_text());c['datasets']=[d for d in c['datasets'] if d['id'] not in ['shot6606','shot6607']]+[{'id':f'shot{s}','label':f'Shot {s} · {date}','path':f'datasets/shot{s}/'} for s,date in [(6606,'2025-02-24'),(6607,'2025-02-25')]];p.write_text(json.dumps(c,indent=2)+'\n')
