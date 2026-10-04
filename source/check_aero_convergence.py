from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1];path=ROOT/'source/analyze_glider.py'
prefix=path.read_text(encoding='utf-8').split('rows=[];last=None;solutions={}')[0]
scope={'__file__':str(path)};exec(compile(prefix.replace('N=14;n=','N=20;n='),str(path),'exec'),scope)
row,_,_=scope['solve'](7.5)
original=json.loads((ROOT/'results/engineering_results.json').read_text())['cruise']
keys=['CL_wing','CD_total','L_D','elevator_trim_deg','static_margin']
result={key:{'N14':float(original[key]),'N20':float(row[key]),'absolute_change':float(row[key]-original[key])} for key in keys}
(ROOT/'results/aero_convergence.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
