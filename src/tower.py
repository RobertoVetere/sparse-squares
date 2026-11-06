"""Exact tower experiments: P_k = prod_i p(x^{m^i}); report dense?, |P_k|, |P_k^2|, exponent log|P^2|/log|P|."""
from flint import fmpz_poly
import math, sys
def stretch(p, m):
    out = [0] * (m * (len(p) - 1) + 1)
    for i, c in enumerate(p): out[m * i] = c
    return fmpz_poly(out)
def tower(p, m, K):
    P = fmpz_poly([1]); rows = []
    for k in range(K):
        P = P * stretch(p, m ** k)
        co = P.coeffs(); sq = (P * P).coeffs()
        t = sum(1 for c in co if c); s = sum(1 for c in sq if c)
        rows.append((k + 1, len(co) - 1, all(co), t, s, math.log(s) / math.log(t)))
    return rows
if __name__ == "__main__":
    import json
    p = json.loads(sys.argv[1]); m = int(sys.argv[2]); K = int(sys.argv[3])
    for r in tower(p, m, K): print("k=%d deg=%d dense=%s |P|=%d |P^2|=%d exp=%.4f" % r)
