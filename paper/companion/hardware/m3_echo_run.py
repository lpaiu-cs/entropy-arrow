"""Phase 3 (B) -- run M3-echo on real IBM hardware: the perspectival control.

On the same scrambled state, contrast a PASSIVE read of the fact site q0 (forward U only -> the
fact has left q0) with a REVERSIBLE Loschmidt echo (forward U, then U^dagger, then read q0 ->
recovered). Recovery is read directly from the R-q0 Bell correlators <ZZ>, <XX> (=+1 each for a
Bell pair, I2~2; ~0 if lost) -- NO randomized measurements, so only 4 circuits per depth
(passive/echo x Z/X). Budget-cheap; the challenge is the echo's 2t gate depth.

Job mode (SamplerV2(mode=backend)) -> Open-Plan compatible. Validate with 0 QPU:
  python m3_echo_run.py --dry-run          # local Aer, Heron-like noise
  python m3_echo_run.py --fake FakeAachen  # full hardware code path, no QPU
  python m3_echo_run.py --backend ibm_marrakesh
"""

from __future__ import annotations
import argparse, json
import numpy as np

try:
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
except Exception as e:
    print("qiskit not available -- hardware artifact; install to run.", e)

from m3_echo_qiskit import m3_echo_circuit, correlator, recovery_MI, NATIVE_BASIS
from m1_ibm_run import pick_line

CFG = dict(N=7, depths=[0, 2, 4, 6, 8, 10], shots=4096, circ_seed=1)


def build_jobs(cfg):
    """(depth, echo, basis, circuit) for the passive/echo x Z/X contrast."""
    jobs = []
    for t in cfg["depths"]:
        for echo in (False, True):
            for basis in ("Z", "X"):
                qc = m3_echo_circuit(cfg["N"], t, cfg["circ_seed"], basis=basis, echo=echo)
                jobs.append((t, echo, basis, qc))
    return jobs


def analyze(counts_map, depths):
    """counts_map[(t,echo,basis)] -> counts. Returns recovery(t) for passive and echo."""
    out = {"depths": depths, "passive": [], "echo": []}
    for t in depths:
        for echo, key in ((False, "passive"), (True, "echo")):
            zz = correlator(counts_map[(t, echo, "Z")])
            xx = correlator(counts_map[(t, echo, "X")])
            out[key].append(recovery_MI(zz, xx))
    out["passive"] = list(map(float, out["passive"]))
    out["echo"] = list(map(float, out["echo"]))
    return out


def run(cfg, backend=None, dry=False):
    jobs = build_jobs(cfg)
    print(f"{len(jobs)} circuits ({len(cfg['depths'])} depths x passive/echo x Z/X), "
          f"{cfg['shots']} shots each")
    counts_map = {}
    if dry:
        from qiskit_aer import AerSimulator
        from qiskit_aer.noise import NoiseModel, depolarizing_error
        nm = NoiseModel()
        nm.add_all_qubit_quantum_error(depolarizing_error(0.003, 2), ["cz"])
        nm.add_all_qubit_quantum_error(depolarizing_error(0.0003, 1), ["sx", "x", "u"])
        sim = AerSimulator(noise_model=nm)
        for (t, e, b, qc) in jobs:
            counts_map[(t, e, b)] = sim.run(qc.decompose(), shots=cfg["shots"]).result().get_counts()
        return counts_map

    from qiskit_ibm_runtime import SamplerV2
    line = pick_line(backend, cfg["N"] + 1)
    is_local = "fake" in backend.name.lower()
    layout = list(line[1:]) + [line[0]]              # R=line[0] adjacent to q0=line[1]; chain on line[1:]
    print(f"backend {backend.name}; line {line} (R@{line[0]} next to q0@{line[1]})  "
          f"({'LOCAL/fake' if is_local else 'HARDWARE'})", flush=True)
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1, initial_layout=layout)
    isa = [pm.run(qc) for (_, _, _, qc) in jobs]
    sampler = SamplerV2(mode=backend)                    # job mode: Open-Plan compatible
    try:
        sampler.options.twirling.enable_measure = True
    except Exception:
        pass
    print(f"  submitting job with {len(isa)} circuits ...", flush=True)
    res = sampler.run([(c,) for c in isa], shots=cfg["shots"]).result()
    for (t, e, b, _), r in zip(jobs, res):
        counts_map[(t, e, b)] = r.data.c.get_counts()
    return counts_map


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--fake", default=None)
    ap.add_argument("--backend", default=None)
    ap.add_argument("--instance", default=None)
    ap.add_argument("--shots", type=int, default=CFG["shots"])
    ap.add_argument("--out", default="m3_result_hw.json")
    a = ap.parse_args()
    cfg = dict(CFG, shots=a.shots)

    if a.dry_run:
        print("=== DRY RUN (local Aer, Heron-like noise) ===")
        cm = run(cfg, dry=True)
    elif a.fake:
        print(f"=== FAKE VALIDATION ({a.fake}, no QPU) ===")
        import qiskit_ibm_runtime.fake_provider as fp
        cm = run(cfg, backend=getattr(fp, a.fake)())
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        s = QiskitRuntimeService(instance=a.instance) if a.instance else QiskitRuntimeService()
        backend = s.backend(a.backend) if a.backend else \
            s.least_busy(operational=True, simulator=False, min_num_qubits=cfg["N"] + 1)
        cm = run(cfg, backend=backend)

    res = analyze(cm, cfg["depths"])
    print("\n--- M3-echo results (recovery I2(R:q0), 0..2 bits) ---")
    print(f"{'depth':>6} {'passive':>8} {'echo':>6} {'gap':>6}")
    for i, t in enumerate(cfg["depths"]):
        p, e = res["passive"][i], res["echo"][i]
        print(f"{t:>6} {p:>8.2f} {e:>6.2f} {e-p:>6.2f}")
    print("\n  passive -> 0 (fact leaves q0); echo stays high (reversible recovery) until 2t noise "
          "wins.\n  The gap is the perspectival signal: same record, recovered iff the observer is "
          "reversible.")
    json.dump(res, open(a.out, "w"), indent=2)
    print(f"saved {a.out}")


if __name__ == "__main__":
    main()
