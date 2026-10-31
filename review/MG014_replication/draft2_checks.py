"""Static MG5 Draft 2 review checks; no MG5 trajectory or event list generated.

Sparse-record checks exercise the existing MG4 detector, not an independent
implementation. Their purpose is to check whether omitted nodes can matter.
"""
from pathlib import Path
import json, sys
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'inputs' / 'src'))
import measure_graph as mg
import graphs as gr


def sparse_check(label, g, s, threshold):
    full = mg.particles(g, s, threshold)
    keep = np.flatnonzero(s < -1.0)  # preserves original row ordering / tie break
    small = dict(pos=g['pos'][keep], N=len(keep), Lx=g['Lx'], Ly=g['Ly'])
    reduced = mg.particles(small, s[keep], threshold) if len(keep) else []
    assert len(full) == len(reduced)
    largest = 0.0
    for a, b in zip(full, reduced):
        assert a['node'] == int(keep[b['node']])
        assert a['depth'] == b['depth']
        err = float(np.max(np.abs(gr.min_image(a['centre'] - b['centre'], g['Lx'], g['Ly']))))
        largest = max(largest, err)
        assert err <= 1e-12
    return dict(case=label, threshold=threshold, full_nodes=len(s), saved_nodes=len(keep),
                particle_count=len(full), max_centre_coordinate_difference=largest)


def main():
    out = dict(commit='a1c4b14c81cd5f8c72956b4ac7da9cc05153d491',
               proposal_blob='470d5621e12fb689e1444d83d5011266afaf03d0')
    n = 18901
    rows = []
    for part, nu, horizon, runs, tick in [
            ('A slow', .002, 500, 4, 1), ('A medium', .02, 500, 4, 1),
            ('A fast', .2, 500, 4, .2), ('B fast', .2, 300, 4, .2)]:
        births = n * nu * horizon
        events = 2 * births
        # Illustrative uncompressed binary schema: float64 time, uint8 type,
        # uint64 ID per event, plus two float64 coordinates per birth.
        byte_count = 17 * events + 16 * births
        rows.append(dict(group=part, runs=runs, expected_births_per_run=births,
                         expected_events_per_run=events,
                         expected_events_per_tick=2*n*nu*tick,
                         expected_event_bytes_per_run=byte_count))
    out['event_storage'] = rows
    out['total_expected_events'] = sum(r['runs']*r['expected_events_per_run'] for r in rows)
    out['total_expected_event_bytes_schema25'] = sum(r['runs']*r['expected_event_bytes_per_run'] for r in rows)
    # uint32 IDs, no alignment padding: 13 bytes/event plus 16/birth.
    out['total_expected_event_bytes_schema21'] = 21*out['total_expected_events']
    assert round(out['total_expected_events']) == 25856568
    assert round(out['total_expected_event_bytes_schema25']) == 646414200
    out['aligned_tick_step_combinations'] = []
    for tick in [1., .2, .1, .05]:
        for dt in [.05, .025]:
            steps = round(tick/dt)
            assert abs(tick/dt-steps) < 1e-12
            assert abs(10/tick-round(10/tick)) < 1e-12
            out['aligned_tick_step_combinations'].append(dict(tick=tick, dt=dt, steps=steps))

    # Scalar batch example: one old node dies, one birth survives, one is born
    # and dies in the tick. Both bookkeeping conventions work; mixing fails.
    before, removed_old, surviving_birth, ephemeral = 10., 2., 3., 7.
    after = before-removed_old+surviving_birth
    assert after == before-(removed_old+ephemeral)+(surviving_birth+ephemeral)
    out['batch_jump_example'] = dict(before=before, after=after,
        effective=dict(R=removed_old, B=surviving_birth),
        all_events=dict(R=removed_old+ephemeral, B=surviving_birth+ephemeral),
        mixed_convention_residual=after-(before-(removed_old+ephemeral)+surviving_birth))

    checks = []
    import run_mg4c as rm
    dec = json.loads((ROOT/'inputs/results/MG4/part0/decisions.json').read_text())
    for slot, pl in [(0,0), (0,1), (1,0), (1,1)]:
        job = rm.job(28, 'D', slot, pl)
        g, _, _ = rm.graph_for(job, dec)
        state = np.load(ROOT/f"inputs/results/MG4c/fields/{job['key']}.npz")['final']
        s = mg.smoother(g) @ state.reshape(3,-1)[1]
        checks.append(sparse_check(job['key'], g, s, -4.457582))
    # Synthetic already-smoothed fields: periodic wrap, exact equal-depth tie,
    # multiple particles, and no candidate. No ODE integrations involved.
    g = dict(pos=np.array([[.05,.1],[7.95,.1],[.1,.2],[4.,4.],[4.1,4.],
                          [2.,2.],[6.,6.],[7.,3.]]), N=8, Lx=8., Ly=8.)
    for label, s in [('wrap_tie_multiple', [-8,-8,-6,-7,-6,0,-2,1]),
                     ('no_candidate', [0,0,-2,0,0,0,0,0]),
                     ('empty_sparse', [0,0,0,0,0,0,0,0])]:
        for threshold in [-4.457582, -5.117354, -1.]:
            checks.append(sparse_check(label, g, np.array(s,dtype=float), threshold))
    out['sparse_record_checks'] = checks
    out['sparse_record_check_scope'] = ('Existing MG4 detector; same float64 smoothed values; '
        'original relative row order preserved; centroids agree within 1e-12, not claimed bitwise identical. '
        'Does not independently verify smoothing, ROS2, interpolation or any MG5 outcome.')
    (ROOT/'draft2_checks.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
