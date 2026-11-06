"""MG10 static design arithmetic only: reads a committed reference; no dynamics."""
from pathlib import Path
import hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[1]
def calculate():
    p=ROOT/"results/MG003/ref_profiles.npz"
    key="h0.2_b29.0"
    with np.load(p) as z: phi=np.array(z[key][1:],float)
    rms=np.sqrt(np.mean(phi**2,axis=1))
    steady=np.array([140/9,140/9,29/(140/9)])
    grid=[1,10,100,1000,10000,100000,1000000]
    rows=[]
    for C in grid:
        depletion=float(-np.expm1(-1/C))
        rows.append(dict(capacity=C,A_fraction_consumed_at_Tref=depletion,
          homogeneous_quasistatic_G_discrepancy=depletion*float(steady[0]),
          G_error_normalized_by_reference_pattern=depletion*float(steady[0]/rms[0]),
          initial_A_per_Omega=C*14*300,
          initial_B_per_Omega=C*29*float(steady[1])*300))
    return dict(kind="STATIC_ARITHMETIC_NOT_A_SIMULATION",reference_key=key,
      reference_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
      reference_perturbation_rms=rms.tolist(),steady_concentrations=steady.tolist(),
      T_ref=300,rows=rows,
      approximate_capacity_for_one_percent_G=float(steady[0]/(.01*rms[0])),
      caveat="A(t)/A0=exp(-t/(C*T_ref)) is exact for deterministic P17 A depletion. G discrepancy uses the homogeneous quasistatic response G=a*r_A/(q-g*p), not a solved finite-bath trajectory or a rigorous soliton-error bound.",
      original_grid_primary_runs=12,extended_grid_primary_runs=16,
      original_grid_with_three_Euler_steps_and_one_ODE_reference=48,
      extended_grid_with_three_Euler_steps_and_one_ODE_reference=64,
      conservative_raw_field_bytes=64*301*7*500*8)
if __name__=="__main__":
    print(json.dumps(calculate(),indent=2))
