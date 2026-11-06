"""MG042 all-capacity/full-horizon exact-case check on ONE homogeneous node.

Not one of the N=500 scientific trajectories. Raw results retain every state;
the independent checker verifies availability and accounting without model imports.
"""
import sys,time,json
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import mg10_reservoir as model
from check_mg10 import audit_raw


def run():
    started=time.perf_counter()
    directory=model.ROOT/"results/MG10/qualification/long_availability"
    directory.mkdir(exist_ok=True)
    rows=[]
    for method in model.STEPS:
        for cap in model.CAPACITIES:
            job=dict(id=f"development_one_node_{method}_C{cap}",mode="FINITE",start="H",method=method,capacity=cap,T=300,n=1)
            times,z=model.trajectory(job)
            file=directory/(job["id"]+".npz")
            np.savez_compressed(file,times=times,state=z)
            with np.load(file) as raw: audit=audit_raw(raw["times"],raw["state"],job)
            assert audit["passed"],str((job,audit))
            continuous=np.exp(-times/(cap*300))
            rows.append(dict(job=job,audit=audit,continuous_discretization_error=float(np.max(abs(1-z[:,3,0]-continuous)))))
    result=dict(status="PASS",scope="ONE_NODE_DEVELOPMENT_ONLY",rows=rows,elapsed_worker_seconds=time.perf_counter()-started)
    (directory/"summary.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print({k:result[k] for k in ("status","scope","elapsed_worker_seconds")})


if __name__=="__main__": run()
