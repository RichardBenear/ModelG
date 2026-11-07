"""MG052 synthetic field/supply measurements; no trajectory generation."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import json
import numpy as np
from mg10d1 import fixture, measurements, ROOT


def main():
    case=fixture(2,0)
    data=np.zeros((1,1,22,500),dtype=np.int64)
    data[0,0,:7]=case["initial"]
    checks=[]
    for label,value in (("flat_original",15556),("flat_depleted_below_old_threshold",1000)):
        data[0,0,1]=value
        result=measurements(data,case)
        assert np.all(result["field_rms"]==0)
        assert abs(result["field_mean"][0,0,1]-value/1000)<1e-12
        checks.append(dict(name=label,passed=True))
    data[0,0,1]=np.tile([9000,11000],250)
    data[0,0,3]=np.tile([0,int(case["a0"][0])],250)
    data[0,0,4]=int(case["b0"][0])//2
    data[0,0,7:14]=np.arange(7)[:,None]
    data[0,0,20:22]=np.array([3,7])[:,None]
    result=measurements(data,case)
    assert result["field_mean"][0,0,1]==10 and result["field_rms"][0,0,1]==1
    assert result["availability_mean"][0,0,0]==.5
    assert result["availability_sd"][0,0,0]==.5
    assert result["availability_min"][0,0,0]==0 and result["availability_max"][0,0,0]==1
    expected_b=(int(case["b0"][0])//2)/case["b0"][0]
    assert abs(result["availability_mean"][0,0,1]-expected_b)<1e-14
    assert np.array_equal(result["reaction_totals"][0,0],np.arange(7)*500)
    assert np.array_equal(result["external_totals"][0,0],[1500,3500])
    checks.append(dict(name="known_alternating_field_supply_and_channel_totals",passed=True))
    output=dict(status="SYNTHETIC MEASUREMENTS PASS",checks=checks,
                scope="Field mean/RMS and stock/channel summaries only; no particle localization or identity qualification")
    (ROOT/"results/MG10D1/development/measurement_check.json").write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(output,indent=2))


if __name__=="__main__":main()
