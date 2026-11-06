"""Replayable actual owned-process stop checks; output directory must be new."""
import sys,json,time
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from run_mg10 import execute


def run(base):
    base=Path(base);base.mkdir(parents=True,exist_ok=False)
    began=time.perf_counter()
    job=dict(id="development_stop",mode="FINITE",start="H",method="E2",capacity=1,T=10,n=8)
    rows=[]
    for name,ceiling,timeout in [("timeout",30,.05),("budget",.05,30)]:
        folder=base/name
        ledger=execute([job],folder,ceiling,timeout,1)
        assert len(ledger["records"])==1
        status=ledger["records"][0]["status"]
        assert status==("TIMEOUT" if name=="timeout" else "BUDGET-SHORT")
        resumed=execute([job],folder,ceiling,timeout,1)
        assert resumed["charged_worker_seconds"]==ledger["charged_worker_seconds"]
        assert len(resumed["records"])==1
        rows.append(dict(check=name,status=status,cost=ledger["charged_worker_seconds"]))
    try:
        with patch("run_mg10.time.sleep",side_effect=KeyboardInterrupt):
            execute([job],base/"interruption",30,30,1)
    except KeyboardInterrupt:
        ledger=json.loads((base/"interruption/ledger.json").read_text())
        assert ledger["records"][0]["status"]=="INTERRUPTED"
        rows.append(dict(check="interruption",status="INTERRUPTED",cost=ledger["charged_worker_seconds"]))
    else: raise AssertionError("No interruption")
    result=dict(status="PASS",rows=rows,elapsed_worker_seconds=time.perf_counter()-began)
    (base/"summary.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(result)


if __name__=="__main__": run(sys.argv[1])
