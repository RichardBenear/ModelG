"""Arguments: path to ModelG checkout with src and MG4 Part 0, path to results/MG4c."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(sys.argv[1]) / "src"))
import evolve
import mg6
import run_mg4c

archive = Path(sys.argv[2])
decisions = json.loads(Path(run_mg4c.DEC).read_text())
tables = json.loads((archive / "evaluation.json").read_text())["tables"]
rows = []
for slot in range(2):
    graph, _, _ = run_mg4c.graph_for(dict(kind="D", slot=slot, pl=0), decisions)
    lx, ly = graph["Lx"], graph["Ly"]
    pos = np.asarray(graph["pos"])
    mesh = evolve.Mesh(pos, lx, ly)
    grid = mg6.Grid(lx, ly)
    for pl in range(2):
        for b in [27.5, 27.75, 28., 28.25, 28.5]:
            key = f"b{b:.6f}_D{slot}_p{pl}"
            with np.load(archive / "fields" / (key + ".npz")) as z:
                field = z["final"]
            displacement = pos - pos[field[1].argmin()]
            periods = np.array([lx, ly])
            radius = np.linalg.norm(displacement-periods*np.rint(displacement/periods),axis=1)
            nodal = np.expm1(2*field[0, radius < 1.5]).sum()
            smooth = evolve.smoothed_X(graph, field)
            record = dict(key=key, b=b, node_P=nodal, max_G=field[0].max(),
                          node_M=np.maximum(-field[0],0).mean()*grid.A)
            interpolated = grid.field(mesh, field)
            for label, threshold in [("common", tables["b28.000000"]["theta"]),
                                     ("by_b", tables[f"b{b:.6f}"]["theta"])]:
                particles = evolve.detect_full(pos, smooth, threshold, lx, ly)
                assert len(particles) == 1
                weights = grid.disc_weights(particles[0]["centre"], 1.5)
                record[label+"_P"] = len(pos)/grid.A * np.dot(weights,np.expm1(2*interpolated))
            rows.append(record)
print(json.dumps(rows, indent=2))
