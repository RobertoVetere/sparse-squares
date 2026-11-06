#!/usr/bin/env python3
"""Independent check of the overlapping tower (does not import src/).
P_K(x) = prod_{j<K} p(L * x^(m^j)),  p = certified seed (midpoint, 100 digits), L = scaling prime.
Checks for small K: (1) P_K dense, (2) |P_K^2| computed directly in high precision,
(3) digit-sum bound |{sum s_j m^j : s_j in S}| by plain set enumeration, (4) growth ratios of the bound."""
import json, sys, math
import mpmath as mp
mp.mp.dps = 120
d = json.load(open(sys.argv[1])); n, m, low = d["n"], d["m"], d["low"]; L = int(sys.argv[2]); Kmax = int(sys.argv[3])
# re-derive high-precision coefficients by Newton on the full (non-reduced) equations, independently
h = n // 2; Z = [j for j in range(2, n + 1) if j not in low]
def full(u):
    a = [mp.mpf(1)] + list(u)
    for k in range(h + 1, n + 1): a.append(d["eps"] * d["sg"] ** (n - k) * a[n - k])
    return a
def F(u):
    a = full(u); return mp.matrix([mp.fsum(a[i] * a[j - i] for i in range(j + 1)) for j in Z])
u = [mp.mpf(x) for x in d["a"][1:h + 1]]
for _ in range(80):
    f0 = F(u); J = mp.matrix(len(Z), h); e = mp.mpf(10) ** -60
    for l in range(h):
        up = list(u); up[l] += e; J[:, l] = (F(up) - f0) / e
    du = mp.lu_solve(J, -f0); u = [u[i] + du[i] for i in range(h)]
    if mp.norm(du) < mp.mpf(10) ** -100: break
a = full(u)
S = [j for j in range(2 * n + 1) if abs(mp.fsum(a[i] * a[j - i] for i in range(max(0, j - n), min(j, n) + 1))) > mp.mpf(10) ** -80]
print("seed n=%d m=%d  |p|=%d  |p^2|=%d  S=%s" % (n, m, n + 1, len(S), S))
sums = {0}; prev = 1
for K in range(1, 7):
    sums = {s + t * m ** (K - 1) for s in sums for t in S}
    print("  digit-sum bound K=%d: %d   ratio to K-1: %.4f" % (K, len(sums), len(sums) / prev)); prev = len(sums)
P = [mp.mpf(1)]
for K in range(1, Kmax + 1):
    st = m ** (K - 1); q = {i * st: a[i] * mp.mpf(L) ** i for i in range(n + 1)}
    R = [mp.mpf(0)] * (len(P) + n * st)
    for x, cx in enumerate(P):
        for y, cy in q.items(): R[x + y] += cx * cy
    P = R
    dense = all(abs(c) > 0 for c in P)
    sq = [mp.mpf(0)] * (2 * len(P) - 1)
    for i, x in enumerate(P):
        for j, y in enumerate(P): sq[i + j] += x * y
    mx = max(abs(v) for v in sq); nz = sum(1 for v in sq if abs(v) > mx * mp.mpf(10) ** -60)
    bound = len({sum(t) for t in [[0]]}) if False else None
    print("K=%d deg=%d terms=%d dense=%s |P^2|=%d  exponent=%.5f  sqrt(deg)=%.1f" % (K, len(P) - 1, len(P), dense, nz, math.log(nz) / math.log(len(P)), math.sqrt(len(P) - 1)), flush=True)
