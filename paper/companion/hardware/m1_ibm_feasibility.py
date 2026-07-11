"""Phase 1-2: device-specific feasibility of M1 on IBM superconducting hardware.

The abstract-brickwork feasibility (paper Fig. 4) used random 2q Cliffords; on IBM each compiles
to ~3 native CZ, tripling the error budget. The hardware-honest choice is a HARDWARE-EFFICIENT
brickwork -- each 2q gate = one native CZ + random single-qubit SU(2) -- a good scrambler at
minimal native cost. This sim rebuilds M1 with that circuit under IBM-realistic noise and the
RENYI-2 mutual information (the quantity randomized measurements actually estimate), to answer:
at IBM error rates, what (N, k, depth) resolves the passive horizon t*_p(k) ~ k/v_B before
decoherence destroys the global record?

Exact density-matrix sim (small N so noise is exact). Noise: 2-qubit depolarizing p2 after each
CZ, single-qubit depolarizing p1 after each SU(2). Readout error + finite shots handled in the
shot-budget section (analytic). A linear chain embeds in IBM heavy-hex with NO routing SWAPs, so
CZ count is not inflated by routing.
"""

import numpy as np

rng = np.random.default_rng(3)
CZ = np.diag([1, 1, 1, -1]).astype(complex)


def haar1():
    z = (rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))) / np.sqrt(2)
    q, r = np.linalg.qr(z); return q * (np.diag(r) / np.abs(np.diag(r)))


def apply_2q(rho, U, i, j, n):
    U4 = U.reshape(2, 2, 2, 2); T = rho.reshape([2] * n + [2] * n)
    T = np.tensordot(U4, T, axes=([2, 3], [i, j])); T = np.moveaxis(T, [0, 1], [i, j])
    T = np.tensordot(T, U4.conj(), axes=([n + i, n + j], [2, 3])); T = np.moveaxis(T, [-2, -1], [n + i, n + j])
    return T.reshape(2 ** n, 2 ** n)


def apply_1q(rho, U, q, n):
    T = rho.reshape([2] * n + [2] * n)
    T = np.tensordot(U, T, axes=([1], [q])); T = np.moveaxis(T, 0, q)
    T = np.tensordot(T, U.conj(), axes=([n + q], [1])); T = np.moveaxis(T, -1, n + q)
    return T.reshape(2 ** n, 2 ** n)


X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]], complex); Z = np.diag([1, -1]).astype(complex)


def depol1(rho, q, n, p):
    if p <= 0: return rho
    return (1 - p) * rho + (p / 3) * (apply_1q(rho, X, q, n) + apply_1q(rho, Y, q, n) + apply_1q(rho, Z, q, n))


def depol2(rho, i, j, n, p):
    if p <= 0: return rho
    acc = rho * (1 - p)
    for P in (X, Y, Z, np.eye(2)):
        for Q in (X, Y, Z, np.eye(2)):
            if P is np.eye(2) and Q is np.eye(2): continue
            acc = acc + (p / 15) * apply_1q(apply_1q(rho, P, i, n), Q, j, n)
    return acc


LET = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
def ptrace(rho, keep, n):
    keep = sorted(keep); T = rho.reshape([2] * n + [2] * n)
    bra = list(LET[:n]); ket = list(LET[n:2 * n])
    for q in range(n):
        if q not in keep: ket[q] = bra[q]
    out = [bra[q] for q in keep] + [ket[q] for q in keep]
    m = 2 ** len(keep)
    return np.einsum(''.join(bra) + ''.join(ket) + '->' + ''.join(out), T, optimize=True).reshape(m, m)


def renyi2(rho):
    return -np.log2(max(np.real(np.trace(rho @ rho)), 1e-15))


def I2(rho, A, B, n):
    return renyi2(ptrace(rho, A, n)) + renyi2(ptrace(rho, B, n)) - renyi2(ptrace(rho, sorted(set(A) | set(B)), n))


def run(Nsys, depth, p2, p1=None, seed=3, ks=(1, 2, 3, 4)):
    global rng; rng = np.random.default_rng(seed)
    p1 = p1 if p1 is not None else p2 / 10.0     # 1q error ~ 0.05% when 2q ~ 0.5%
    n = Nsys + 1; R = Nsys
    psi = np.zeros(2 ** n, complex)
    for b in (0, 1): psi[(b << (n - 1)) | (b << (n - 1 - R))] += 1 / np.sqrt(2)
    rho = np.outer(psi, psi.conj())
    curves = {k: [I2(rho, [R], list(range(k)), n)] for k in ks}
    Itot = [I2(rho, [R], list(range(Nsys)), n)]
    for L in range(1, depth + 1):
        for i in range(L % 2, Nsys - 1, 2):
            rho = apply_1q(rho, haar1(), i, n);     rho = depol1(rho, i, n, p1)
            rho = apply_1q(rho, haar1(), i + 1, n); rho = depol1(rho, i + 1, n, p1)
            rho = apply_2q(rho, CZ, i, i + 1, n);   rho = depol2(rho, i, i + 1, n, p2)
        for k in ks: curves[k].append(I2(rho, [R], list(range(k)), n))
        Itot.append(I2(rho, [R], list(range(Nsys)), n))
    return np.arange(depth + 1), {k: np.array(v) for k, v in curves.items()}, np.array(Itot)


def tstar(ts, y, th=1.0):
    for i in range(len(y) - 1):
        if y[i] >= th > y[i + 1]:
            return ts[i] + (th - y[i]) / (y[i + 1] - y[i]) * (ts[i + 1] - ts[i])
    return np.nan


if __name__ == "__main__":
    import time
    Nsys, depth, ks = 9, 30, (1, 2, 3, 4)
    print("=== M1 on IBM-efficient brickwork (1 CZ/gate), Renyi-2 MI, exact noisy DM ===")
    print("IBM 2q error: Heron ~0.3%, Eagle ~1%. Sweeping p2.")
    OUT = r"C:/Users/lpaiu/AppData/Local/Temp/claude/E--lab-entropy-arrow/30a5eb66-4513-44dd-96c0-3d0d45b958a5/scratchpad"
    store = {}
    for p2 in (0.0, 0.003, 0.005, 0.01):
        t0 = time.time()
        ts, curves, Itot = run(Nsys, depth, p2, ks=ks)
        tp = {k: tstar(ts, curves[k]) for k in ks}
        store[p2] = dict(ts=ts, curves=curves, Itot=Itot, tp=tp)
        # v_B from clean t*_p(k) slope
        kk = np.array(ks); tpv = np.array([tp[k] for k in ks])
        good = np.isfinite(tpv)
        vB = (1 / np.polyfit(kk[good], tpv[good], 1)[0]) if good.sum() >= 2 else np.nan
        print(f"[p2={p2*100:.1f}%] t*_p(k)={{ {', '.join(f'{k}:{tp[k]:.1f}' for k in ks)} }}  "
              f"v_B~{vB:.2f}  I2(R:sys)@end={Itot[-1]:.2f}  ({time.time()-t0:.0f}s)", flush=True)
    np.savez(OUT + "/m1_ibm.npz",
             **{f"tp_{int(p2*1000)}": np.array([store[p2]['tp'][k] for k in ks]) for p2 in store},
             **{f"Itot_{int(p2*1000)}": store[p2]['Itot'] for p2 in store},
             ks=np.array(ks), ts=store[0.0]['ts'])
    print("saved m1_ibm.npz")
