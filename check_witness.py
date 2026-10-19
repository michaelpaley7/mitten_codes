"""Verify a released schedule fault-set witness using the authors' files and Stim."""
import argparse, json
import numpy as np
import stim


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('code_dir'); ap.add_argument('witness')
    ap.add_argument('--cycles', type=int, default=3); ap.add_argument('--fault-cycle', type=int, default=1)
    ap.add_argument('--drop-one', action='store_true')
    a = ap.parse_args()
    assert 0 <= a.fault_cycle < a.cycles
    HX, HZ = np.load(a.code_dir + '/Hx.npy').astype(int), np.load(a.code_dir + '/Hz.npy').astype(int)
    LX, LZ = np.load(a.code_dir + '/Lx.npy').astype(int), np.load(a.code_dir + '/Lz.npy').astype(int)
    S = json.load(open(a.code_dir + '/hook_free_SE_cycle_schedule.json'))
    W = json.load(open(a.witness))
    kind, data_f, anc_f = W['half'], list(W['faults']['data_qubits']), list(W['faults']['ancilla_faults'])
    if a.drop_one:
        (anc_f.pop() if anc_f else data_f.pop())
    n, mx, mz = HX.shape[1], HX.shape[0], HZ.shape[0]
    layers = {k: [[(c, q) for mv in L['moves'] for c, q in mv['cz_gates']] for L in S[f'{k}_layers']['layers']]
              for k in 'XZ'}
    # Each check's gate sequence in the JSON order must reproduce its row of Hx / Hz exactly.
    seq = {k: [[] for _ in range(H.shape[0])] for k, H in (('X', HX), ('Z', HZ))}
    when = {}
    for k in 'XZ':
        for t, L in enumerate(layers[k]):
            for c, q in L:
                seq[k][c].append(q)
                if k == kind:
                    when[(c, len(seq[k][c]))] = t          # layer in which this check's gate number len(...) fires
    for k, H in (('X', HX), ('Z', HZ)):
        assert all(sorted(s) == sorted(np.flatnonzero(H[c]).tolist()) for c, s in enumerate(seq[k])), 'JSON != H'
    for f in anc_f:
        assert seq[kind][f['check_row']][f['after_gate_number']:] == f['spreads_to'], f'spread mismatch: {f}'
    xa = [n + i for i in range(mx)]; za = [n + mx + i for i in range(mz)]
    err, anc = ('X_ERROR', xa) if kind == 'X' else ('Z_ERROR', za)
    basis = 'Z' if kind == 'X' else 'X'                  # memory basis that detects this error type

    def build(with_faults):
        c = stim.Circuit(); c.append('R' if basis == 'Z' else 'RX', range(n))
        nm, last = 0, None
        for cyc in range(a.cycles):
            hit = with_faults and cyc == a.fault_cycle
            c.append('RX', xa); c.append('R', za)
            for k in 'XZ':
                if hit and kind == 'X' and k == 'Z' and data_f:
                    c.append(err, data_f, 1.0)             # X data faults, after the X round
                for t, L in enumerate(layers[k]):
                    c.append('CX', [x for cq in L for x in ((xa[cq[0]], cq[1]) if k == 'X' else (cq[1], za[cq[0]]))])
                    if hit and k == kind:
                        for f in anc_f:
                            if when[(f['check_row'], f['after_gate_number'])] == t:
                                c.append(err, [anc[f['check_row']]], 1.0)
            if hit and kind == 'Z' and data_f:
                c.append(err, data_f, 1.0)                 # Z data faults, after the Z round
            c.append('MX', xa); c.append('M', za)
            rec = list(range(nm + mx, nm + mx + mz)) if basis == 'Z' else list(range(nm, nm + mx)); nm += mx + mz
            for i, r in enumerate(rec):
                c.append('DETECTOR', [stim.target_rec(r - nm)] + ([stim.target_rec(last[i] - nm)] if last else []))
            last = rec
        c.append('M' if basis == 'Z' else 'MX', range(n)); d0 = nm; nm += n
        Hb, Lb = (HZ, LZ) if basis == 'Z' else (HX, LX)
        for i in range(Hb.shape[0]):
            c.append('DETECTOR', [stim.target_rec(d0 + j - nm) for j in np.flatnonzero(Hb[i])] + [stim.target_rec(last[i] - nm)])
        for o in range(Lb.shape[0]):
            c.append('OBSERVABLE_INCLUDE', [stim.target_rec(d0 + j - nm) for j in np.flatnonzero(Lb[o])], o)
        return c

    build(False).detector_error_model()                    # raises if any detector/observable is non-deterministic
    det, obs = build(True).compile_detector_sampler().sample(4, separate_observables=True)
    assert (det == det[0]).all() and (obs == obs[0]).all(), 'faulted circuit is not deterministic'
    nf = len(data_f) + len(anc_f)
    ok = (not det[0].any()) and bool(obs[0].any())
    print(json.dumps(dict(code=S['code'], half=kind, faults=nf, data_faults=len(data_f), ancilla_faults=len(anc_f),
                          cycles=a.cycles, fault_cycle=a.fault_cycle, detectors=int(det.shape[1]),
                          detectors_fired=int(det[0].sum()), logicals_flipped=int(obs[0].sum()),
                          result='PASS: undetectable logical fault set' if ok else 'FAIL')))


if __name__ == '__main__':
    main()
