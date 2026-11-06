"""Beam search over RATIONAL dense seeds built as P = A(x) * (v + u x^k), u/v rational chosen as a rational
root of the quadratic making some coefficient of P^2 vanish. Score = best overlapping-tower exponent."""
import math, sys, itertools
from fractions import Fraction
from math import gcd, isqrt
from flint import fmpz_poly
from digits2 import growth
def sq(p): return [int(c) for c in (fmpz_poly(p) ** 2).coeffs()]
def supp(p): return [i for i, c in enumerate(sq(p)) if c != 0]
def score(p):
    n = len(p) - 1; S = supp(p)
    if len(S) > 2 * n: return (1.0, n + 1)
    return min((math.log(growth(S, m)) / math.log(m), m) for m in range(max(2, (n + 2) // 2), n + 2))
def prim(p):
    g = 0
    for c in p: g = gcd(g, int(c))
    return [int(c) // g for c in p]
def israt_sqrt(D):
    if D < 0: return None
    r = isqrt(D); return r if r * r == D else None
def children(A, kmax_extra=1):
    n = len(A) - 1; A2 = sq(A); out = []
    for k in range(2, n + 2):
        # P^2 coefficient at N: v^2 A2[N] + 2uv A2[N-k] + u^2 A2[N-2k]; root t=u/v of  a + 2b t + c t^2 = 0
        g = lambda i: A2[i] if 0 <= i < len(A2) else 0
        ts = set()
        for N in range(0, len(A2) + 2 * k):
            a, b, c = g(N), g(N - k), g(N - 2 * k)
            if c == 0:
                if b != 0: ts.add(Fraction(-a, 2 * b))
                continue
            r = israt_sqrt(b * b - a * c)
            if r is not None:
                ts.add(Fraction(-b + r, c)); ts.add(Fraction(-b - r, c))
        for t in ts:
            if t == 0: continue
            u, v = t.numerator, t.denominator
            B = [v] + [0] * (k - 1) + [u]
            P = prim([int(x) for x in (fmpz_poly(A) * fmpz_poly(B)).coeffs()])
            if len(P) != n + k + 1 or any(c == 0 for c in P): continue
            out.append(P)
    return out
def seeds():
    base = {"V8": [1, -2, -2, -4, -10, -28, -98, -686, 2401], "A6": [1, 2, -2, 4, -10, 50, 125]}
    # truncated sqrt(1+4x) (integer Catalan-type coefficients)
    s = [1, 2, -2, 4, -10, 28, -84, 264, -858, 2860]
    for L in range(3, 10): base["sqrt%d" % L] = s[:L]
    return base
if __name__ == "__main__":
    beam = list(seeds().values()); seen = set(); width = int(sys.argv[1]); rounds = int(sys.argv[2]); maxdeg = int(sys.argv[3])
    allbest = []
    for r in range(rounds):
        cand = []
        for A in beam:
            for P in children(A):
                key = tuple(P)
                if key in seen or len(P) - 1 > maxdeg: continue
                seen.add(key); e, m = score(P); cand.append((e, m, P))
        cand.sort(key=lambda t: t[0]); allbest += cand[:5]
        beam = [c[2] for c in cand[:width]]
        print("round", r, "candidates", len(cand), "best:", ["%.5f(n=%d,m=%d,|p2|=%d)" % (c[0], len(c[2]) - 1, c[1], len(supp(c[2]))) for c in cand[:4]], flush=True)
    allbest.sort(key=lambda t: t[0])
    print("BEST", "%.5f" % allbest[0][0], "m=%d" % allbest[0][1], allbest[0][2])
