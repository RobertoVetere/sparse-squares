"""Numerical (mp 60 digits) checks for the overlapping tower P_K = prod p(x^{m^i}):
density, actual |P_K^2| vs automaton count, and the set of consecutive-coefficient ratios vs forbidden values."""
import mpmath as mp, json, sys, math
from digits2 import count, growth
from superseed import refine
mp.mp.dps = 60
def level(prev, a, m):
    n = len(a) - 1; deg = (len(prev) - 1) * m + n; c = [mp.mpf(0)] * (deg + 1)
    for q, v in enumerate(prev):
        if v:
            for i, ai in enumerate(a): c[m * q + i] += ai * v
    return c
def run(n, m, low, a0, K, conv_upto=3):
    Z = [j for j in range(2, n + 1) if j not in low]
    a, r = refine(a0, n, 1, 1, Z)
    S = [j for j in range(2*n+1) if abs(mp.fsum(a[i]*a[j-i] for i in range(max(0,j-n), min(j,n)+1))) > mp.mpf(10)**-40]
    print("seed residual", mp.nstr(r,3), "S=", S, "lambda=%.5f exp=%.5f" % (growth(S, m), math.log(growth(S, m))/math.log(m)))
    forb = [-a[r + m] / a[r] for r in range(0, n - m + 1)]
    print("forbidden ratios", [mp.nstr(f, 8) for f in forb])
    P = [mp.mpf(1)]
    for k in range(1, K + 1):
        P = level(P, a, m)  # note: P_k(x) = P_{k-1}(x^m) * p(x)  (same polynomial family)
        mn = min(abs(c) for c in P); ratios = [P[i]/P[i-1] for i in range(1, len(P))]
        dist = min(abs(x - f) / (1 + abs(f)) for x in ratios for f in forb)
        line = "K=%d deg=%d dense_min|c|=%s  min rel dist ratios->forbidden=%s" % (k, len(P)-1, mp.nstr(mn, 4), mp.nstr(dist, 4))
        if k <= conv_upto:
            sq = [mp.mpf(0)] * (2 * len(P) - 1)
            for i, x in enumerate(P):
                for j, y in enumerate(P): sq[i + j] += x * y
            mx = max(abs(v) for v in sq)
            nz = sum(1 for v in sq if abs(v) > mx * mp.mpf(10) ** -35)
            line += "  |P^2|=%d automaton=%d  exp=%.4f" % (nz, count(S, m, k), math.log(nz)/math.log(len(P)))
        print(line, flush=True)
if __name__ == "__main__":
    n, m, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    run(n, m, json.loads(sys.argv[4]), json.loads(sys.argv[5]), K, int(sys.argv[6]) if len(sys.argv) > 6 else 3)
