"""Phase 3 -- run M1 on real IBM hardware and analyze it (submission + analysis pipeline).

Measures the passive memory horizon t*_p(k) ~ k/v_B on an IBM backend and extracts the butterfly
velocity. Efficient design: at each scrambling depth t, measure the (k_max+1)-qubit register
R u [0:k_max] in n_bases random single-qubit bases; the Renyi-2 mutual informations I2(R:[0:k])
for ALL k <= k_max are then MARGINALS of the same data (van Enk-Beenakker randomized-measurement
purities), so the circuit count is only (#depths) x (#bases).

Pipeline:
  1. connect + pick a low-error qubit LINE on the backend (no routing SWAPs for a chain).
  2. build M1 circuits (Bell + hardware-efficient scramble + random local bases + measure).
  3. submit with SamplerV2 in a Session, measurement twirling (TREX) on.
  4. estimate I2(R:[0:k])(t) with bootstrap error bars -> t*_p(k) -> v_B = 1/slope.
  5. (optional) OTOC light-cone cross-check of v_B (note: 2t depth, coherence-limited).

Run `--dry-run` first (local AerSimulator with a noise model) to validate the whole pipeline
before spending hardware time. Requires: qiskit>=1.0, qiskit-ibm-runtime>=0.20, qiskit-aer.

  python m1_ibm_run.py --dry-run
  python m1_ibm_run.py --backend ibm_torino --shots 1024 --n-bases 200
"""

from __future__ import annotations
import argparse, json, time
import numpy as np

try:
    from qiskit import QuantumCircuit
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
except Exception as e:
    QuantumCircuit = None
    print("qiskit not available -- this is the Phase-3 artifact; install to run.", e)

# reuse the validated circuit + estimator pieces
from m1_qiskit import scrambler, _rand_su2_angles

CFG = dict(N=11, kmax=4, depths=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18], n_bases=200, shots=512,
           circ_seed=1, meas_seed0=7000)

# budget-friendly preset for the IBM free Open Plan (a first real horizon; noisier MI).
SMOKE = dict(N=8, kmax=3, depths=[0, 2, 4, 6, 8, 10, 12], n_bases=80, shots=256,
             circ_seed=1, meas_seed0=7000)


# ----------------------------------------------------------------- backend + qubit line
def pick_line(backend, length):
    """Greedy lowest-2q-error path of `length` qubits in the backend coupling graph."""
    cmap = backend.coupling_map
    target = backend.target
    def err(a, b):
        try:
            return target["cz"][(a, b)].error
        except Exception:
            try:
                return target["ecr"][(a, b)].error
            except Exception:
                return 1e-2
    adj = {}
    for a, b in cmap.get_edges():
        adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    best = None
    for start in adj:
        path = [start]; used = {start}
        while len(path) < length:
            cur = path[-1]
            nxt = sorted((n for n in adj.get(cur, ()) if n not in used), key=lambda n: err(cur, n))
            if not nxt:
                break
            path.append(nxt[0]); used.add(nxt[0])
        if len(path) == length:
            cost = sum(err(path[i], path[i + 1]) for i in range(length - 1))
            if best is None or cost < best[0]:
                best = (cost, list(path))
    if best is None:
        raise RuntimeError(f"no line of length {length} found on {backend.name}")
    return best[1]


# ----------------------------------------------------------------- circuits
def m1_scan_circuit(N, depth, kmax, circ_seed, meas_seed):
    """Bell(R,q0) + scramble + random local bases on {q0..q_{kmax-1}, R} + measure those."""
    qc = QuantumCircuit(N + 1, kmax + 1)
    R = N
    qc.h(R); qc.cx(R, 0)
    qc.compose(scrambler(N, depth, circ_seed), qubits=range(N), inplace=True)
    probe = list(range(kmax)) + [R]                    # cbit 0..kmax-1 = window, cbit kmax = R
    rng = np.random.default_rng(meas_seed)
    for q in probe:
        qc.u(*_rand_su2_angles(rng), q)
    for c, q in enumerate(probe):
        qc.measure(q, c)
    return qc


# ----------------------------------------------------------------- RM Renyi-2 (marginal)
def _hamming_weight_matrix(nq):
    """W[a,b] = (-2)^(-Hamming(a,b)) over the 2^nq marginal words (nq small)."""
    d = 1 << nq
    a = np.arange(d)
    xor = a[:, None] ^ a[None, :]
    ham = np.array([[bin(int(v)).count("1") for v in row] for row in xor])
    return (-2.0) ** (-ham)


_WCACHE = {}
def _marg_purity(bitlists, positions):
    """Unbiased same-basis purity of the marginal on `positions` (cbit indices), averaged over
    bases. Histogram form (O(4^nq) per basis, exact): with W the Hamming weight matrix,
    Tr(rho^2) = 2^nq/(M(M-1)) [ h.W.h - M ] (the '-M' removes the s=s' self-pairs, W[a,a]=1)."""
    nq = len(positions); d = 1 << nq
    if nq not in _WCACHE:
        _WCACHE[nq] = _hamming_weight_matrix(nq)
    W = _WCACHE[nq]
    est = 0.0; nb = 0
    for arr in bitlists:
        arr = np.asarray(arr)
        if len(arr) < 2:
            continue
        vals = np.zeros(len(arr), dtype=np.int64)
        for j, p in enumerate(positions):
            vals |= ((arr >> p) & 1) << j
        h = np.bincount(vals, minlength=d).astype(float)
        M = h.sum()
        est += (2 ** nq) * (h @ W @ h - M) / (M * (M - 1)); nb += 1
    return est / max(nb, 1)


def I2_of_k(bitlists, kmax, k):
    """I2(R:[0:k]) = S2(R)+S2([0:k])-S2(R u [0:k]) from marginals. cbits: 0..kmax-1 window, kmax=R."""
    S2 = lambda pos: -np.log2(max(_marg_purity(bitlists, pos), 1e-12))
    win = list(range(k)); R = [kmax]
    return S2(R) + S2(win) - S2(win + R)


def tstar(depths, I2curve, th=1.0):
    d = np.asarray(depths, float)
    for i in range(len(d) - 1):
        if I2curve[i] >= th > I2curve[i + 1]:
            return d[i] + (th - I2curve[i]) / (I2curve[i + 1] - I2curve[i]) * (d[i + 1] - d[i])
    return np.nan


def _fit_vB(tp, kmin=2):
    """v_B = 1/slope of t*_p(k) vs k, fitting k>=kmin (small-k has threshold curvature)."""
    kk = np.array([k for k in tp if k >= kmin and np.isfinite(tp[k])], float)
    tv = np.array([tp[k] for k in kk])
    if len(kk) < 2:                       # fall back to all finite k if too few
        kk = np.array([k for k in tp if np.isfinite(tp[k])], float)
        tv = np.array([tp[k] for k in kk])
    if len(kk) < 2:
        return np.nan, 0
    s = np.polyfit(kk, tv, 1)[0]
    return (1 / s if s > 0 else np.nan), len(kk)


def analyze(bitlists_by_depth, depths, kmax, n_boot=200, kmin=2, rng=None):
    """I2(R:[0:k])(t) for k=1..kmax, t*_p(k) with bootstrap error, and v_B = 1/slope of t*_p(k)
    (fit over k>=kmin to avoid the small-k threshold curvature)."""
    rng = rng or np.random.default_rng(0)
    ks = list(range(1, kmax + 1))
    I2 = {k: np.array([I2_of_k(bitlists_by_depth[t], kmax, k) for t in depths]) for k in ks}
    tp = {k: tstar(depths, I2[k]) for k in ks}
    vB, npts = _fit_vB(tp, kmin)
    vB_bs = []
    for _ in range(n_boot):
        tp_bs = {}
        for k in ks:
            curve = []
            for t in depths:
                bl = bitlists_by_depth[t]
                idx = rng.integers(0, len(bl), len(bl))
                curve.append(I2_of_k([bl[i] for i in idx], kmax, k))
            tp_bs[k] = tstar(depths, curve)
        v, m = _fit_vB(tp_bs, kmin)
        if np.isfinite(v):
            vB_bs.append(v)
    vB_err = float(np.std(vB_bs)) if len(vB_bs) > 5 else np.nan
    return dict(I2=I2, tp=tp, vB=vB, vB_err=vB_err, n_fit=npts, kmin=kmin)


# ----------------------------------------------------------------- run (hardware or dry-run)
def _counts_to_bitwords(counts, ncbits):
    words = []
    for bs, c in counts.items():
        b = bs.replace(" ", "")
        word = int(b[::-1], 2)                          # cbit 0 = rightmost
        words.extend([word] * c)
    return np.array(words, dtype=np.int64)


def run(cfg, backend=None, dry=False):
    N, kmax = cfg["N"], cfg["kmax"]
    depths, nb, shots = cfg["depths"], cfg["n_bases"], cfg["shots"]
    # build all circuits (depth x basis)
    jobs = []                                            # (depth, basis, circuit)
    for t in depths:
        for b in range(nb):
            jobs.append((t, b, m1_scan_circuit(N, t, kmax, cfg["circ_seed"], cfg["meas_seed0"] + t * nb + b)))
    print(f"{len(jobs)} circuits ({len(depths)} depths x {nb} bases), {shots} shots each")

    if dry:
        from qiskit_aer import AerSimulator
        from qiskit_aer.noise import NoiseModel, depolarizing_error
        nm = NoiseModel()
        nm.add_all_qubit_quantum_error(depolarizing_error(0.003, 2), ["cz"])
        nm.add_all_qubit_quantum_error(depolarizing_error(0.0003, 1), ["sx", "x", "u"])
        sim = AerSimulator(noise_model=nm)
        results = {}
        for (t, b, qc) in jobs:
            counts = sim.run(qc.decompose(), shots=shots).result().get_counts()
            results.setdefault(t, []).append(_counts_to_bitwords(counts, kmax + 1))
        return results

    # ---- hardware path ----
    # JOB mode (SamplerV2(mode=backend)): works on the free Open Plan (Session mode is premium-only
    # and 400s on Open with "not authorized to run a session"). Batch mode also works if available.
    from qiskit_ibm_runtime import SamplerV2
    is_local = "fake" in backend.name.lower() \
        or backend.__class__.__module__.startswith("qiskit_ibm_runtime.fake_provider")
    line = pick_line(backend, N + 1)
    # logical order is [system 0..N-1, R=N]; the circuit does cx(R, q0) and the scrambler is a
    # chain on 0..N-1, so map the system chain to line[1:] and the reference R to line[0], which
    # is adjacent to q0=line[1] -- no routing SWAP for the Bell prep or the chain.
    layout = list(line[1:]) + [line[0]]
    print(f"backend {backend.name}; line {line} (R@{line[0]} next to q0@{line[1]})  "
          f"({'LOCAL/fake' if is_local else 'HARDWARE'})", flush=True)
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1, initial_layout=layout)
    isa = [pm.run(qc) for (_, _, qc) in jobs]

    sampler = SamplerV2(mode=backend)                      # job mode: Open-Plan compatible
    try:
        sampler.options.twirling.enable_measure = True     # TREX measurement-error mitigation
    except Exception:
        pass
    out, B = [], 300                                       # <=300 PUBs per job
    for s in range(0, len(isa), B):
        print(f"  submitting job with {len(isa[s:s+B])} circuits ...", flush=True)
        res = sampler.run([(c,) for c in isa[s:s + B]], shots=shots).result()
        out.extend(res)

    results = {}
    for (t, b, _), r in zip(jobs, out):
        counts = r.data.c.get_counts()                     # robust across runtime versions
        results.setdefault(t, []).append(_counts_to_bitwords(counts, kmax + 1))
    return results


# ----------------------------------------------------------------- OTOC cross-check (optional)
def otoc_circuit(N, depth, x, circ_seed):
    """Interferometric OTOC of edge V=X_0 and W=Z_x: ancilla-controlled protocol (2t depth).
    Returns a circuit whose <X_anc> gives Re OTOC(x,t). v_B = wavefront slope where OTOC ~ 1/2."""
    qc = QuantumCircuit(N + 1, 1); anc = N
    U = scrambler(N, depth, circ_seed)
    qc.h(anc)
    qc.cx(anc, 0)                                          # controlled V = X_0
    qc.compose(U, range(N), inplace=True)
    qc.z(x)                                                # W = Z_x
    qc.compose(U.inverse(), range(N), inplace=True)
    qc.cx(anc, 0)
    qc.h(anc); qc.measure(anc, 0)
    return qc


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--fake", default=None, help="validate the HARDWARE path locally on a fake "
                    "backend (e.g. FakeAachen); spends no QPU time")
    ap.add_argument("--smoke", action="store_true", help="small budget-friendly preset (free plan)")
    ap.add_argument("--backend", default=None)
    ap.add_argument("--instance", default=None, help="IBM Cloud instance CRN/name (new platform)")
    ap.add_argument("--shots", type=int, default=None)
    ap.add_argument("--n-bases", type=int, default=None)
    ap.add_argument("--out", default="m1_result.json")
    a = ap.parse_args()
    base = SMOKE if a.smoke else CFG
    cfg = dict(base, shots=a.shots or base["shots"], n_bases=a.n_bases or base["n_bases"])
    ncirc = len(cfg["depths"]) * cfg["n_bases"]
    print(f"config: N={cfg['N']} kmax={cfg['kmax']} depths={cfg['depths']} "
          f"n_bases={cfg['n_bases']} shots={cfg['shots']}")
    print(f"  -> {ncirc} circuits x {cfg['shots']} shots = {ncirc*cfg['shots']:,} total shots"
          f"  (free Open Plan is minutes/month of QPU time -- start with --smoke)")

    if a.dry_run:
        print("=== DRY RUN (local AerSimulator, IBM Heron-like noise) ===")
        results = run(cfg, dry=True)
    elif a.fake:
        print(f"=== FAKE-BACKEND VALIDATION ({a.fake}, local, no QPU time) ===")
        import qiskit_ibm_runtime.fake_provider as fp
        backend = getattr(fp, a.fake)()
        results = run(cfg, backend=backend, dry=False)
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        service = QiskitRuntimeService(instance=a.instance) if a.instance else QiskitRuntimeService()
        backend = (service.backend(a.backend) if a.backend
                   else service.least_busy(operational=True, simulator=False,
                                           min_num_qubits=cfg["N"] + 1))
        results = run(cfg, backend=backend, dry=False)

    res = analyze(results, cfg["depths"], cfg["kmax"])
    print("\n--- results ---")
    for k in range(1, cfg["kmax"] + 1):
        print(f"  t*_p(k={k}) = {res['tp'][k]:.1f} CZ   I2(R:[0:{k}])(t) = {np.round(res['I2'][k],2)}")
    print(f"\n  v_B = {res['vB']:.3f} +/- {res['vB_err']:.3f}  (slope fit over {res['n_fit']} "
          f"points, k>={res['kmin']}; cross-check vs OTOC wavefront)")
    json.dump({"tp": {k: float(v) for k, v in res["tp"].items()},
               "vB": float(res["vB"]), "vB_err": float(res["vB_err"]),
               "I2": {k: res["I2"][k].tolist() for k in res["I2"]},
               "depths": cfg["depths"], "cfg": {k: cfg[k] for k in ("N", "kmax", "n_bases", "shots")}},
              open(a.out, "w"), indent=2)
    print(f"saved {a.out}")


if __name__ == "__main__" and QuantumCircuit is not None:
    main()
