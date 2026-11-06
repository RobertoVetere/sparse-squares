"""Super-seed analysis with high precision.
p: real seed (refined to high precision by Newton on the symmetric system), m: overlap base, K: levels.
P_K = prod_{i<K} p(x^{m^i}).  supp(P_K^2) is bounded by digit sums of S=supp(p^2) (exact combinatorics).
Then the non-overlapping tower of P_K with base B=deg P_K+1 is dense and has growth lambda'(S',B)."""
import mpmath as mp, numpy as np, math, json, sys
from digits2 import build
from digits import growth as growth01
mp.mp.dps = 80
def refine(a0, n, eps, sg, Z):
    h = n // 2
    u = [mp.mpf(x) for x in a0[1:h+1]]
    def full(u):
        a = [mp.mpf(1)] + list(u) + [None] * (n - h)
        for k in range(h + 1, n + 1): a[k] = eps * sg ** (n - k) * a[n - k]
        return a
    def F(u):
        a = full(u); return [mp.fsum(a[i] * a[j - i] for i in range(j + 1)) for j in Z]
    for _ in range(60):
        Fv = mp.matrix(F(u)); Jm = mp.matrix(len(Z), h); e = mp.mpf(10) ** (-40)
        for l in range(h):
            up = list(u); up[l] += e; Fp = F(up)
            for r in range(len(Z)): Jm[r, l] = (Fp[r] - Fv[r]) / e
        d = mp.lu_solve(Jm, -Fv); u = [u[i] + d[i] for i in range(h)]
        if mp.norm(d) < mp.mpf(10) ** (-70): break
    return full(u), max(abs(x) for x in F(u))
def digit_support(S, m, K):
    sums = {0}
    for i in range(K): sums = {s + t * m ** i for s in sums for t in S}
    return sums
def analyse(a, n, m, K):
    S = sorted(j for j in range(2 * n + 1) if abs(mp.fsum(a[i] * a[j - i] for i in range(max(0, j - n), min(j, n) + 1))) > mp.mpf(10) ** (-50))
    # build P_K coefficients in high precision
    P = [mp.mpf(1)]
    for i in range(K):
        st = m ** i; q = [mp.mpf(0)] * (n * st + 1)
        for j, c in enumerate(a): q[j * st] = c
        R = [mp.mpf(0)] * (len(P) + len(q) - 1)
        for x, cx in enumerate(P):
            if cx == 0: continue
            for y, cy in enumerate(q):
                if cy != 0: R[x + y] += cx * cy
        P = R
    deg = len(P) - 1; minabs = min(abs(c) for c in P); maxabs = max(abs(c) for c in P)
    Sp = digit_support(S, m, K)  # superset of supp(P_K^2)
    B = deg + 1
    lam = growth01(Sp, B)
    return dict(n=n, m=m, K=K, S=S, deg=deg, terms=deg + 1, dense_minabs=float(minabs), maxabs=float(maxabs),
                sq_support_bound=len(Sp), exp_superseed=math.log(lam) / math.log(B), lam=lam,
                exp_naive=math.log(len(Sp)) / math.log(deg + 1))
if __name__ == "__main__":
    n, m, Kmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); low = json.loads(sys.argv[4]); a0 = json.loads(sys.argv[5])
    Z = [j for j in range(2, n + 1) if j not in low]
    a, r = refine(a0, n, 1, 1, Z); print("residual", mp.nstr(r, 5))
    for K in range(1, Kmax + 1):
        d = analyse(a, n, m, K); print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in d.items() if k != "S"}), flush=True)
