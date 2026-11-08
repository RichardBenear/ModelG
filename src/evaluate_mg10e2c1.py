"""MG076 same frozen E2 measurements, separate continuation report/namespace."""
from pathlib import Path
from evaluate_mg10e2 import main as original_evaluate
from mg10d1_pool import atomic_json
import json
BASE=Path(__file__).resolve().parents[1]/'results/MG10E2C1'
def main(base=BASE):
    base=Path(base);original_evaluate(base)
    p=base/'MG10E2_results.md';text=p.read_text(encoding='utf-8').replace('# MG10-E2 results','# MG10-E2-C1 continuation results',1)
    text+='\nThis is the separately registered continuation epoch. Original E2 remains incomplete as recorded;32 prepared parents and1 completed branch are reused,7 interrupted paths require exact-prefix seed replay,120 original unstarted paths remain in the same fixed128-branch matrix. No original report/data overwritten. Aggregate accounting includes the first attempt; no physical identity/quantum-memory or full-recovery claim.\n'
    p.rename(base/'MG10E2C1_results.md');(base/'MG10E2C1_results.md').write_text(text,encoding='utf-8')
    e=json.loads((base/'evaluation.json').read_text());e['epoch']='MG076 E2-C1';e['original_record']='E2 remains INCOMPLETE RESTORATION ESTIMATION';atomic_json(base/'evaluation.json',e)
if __name__=='__main__':main()
