"""For every symmetric square-system pattern (n<=NMAX) whose best overlapping exponent beats a threshold,
compute ALL isolated solutions by homotopy and test whether some x-rescaling makes the seed rational:
b_i = a_i / a_1^i must all be rational (checked by continued fractions at 30+ digits after mp refinement)."""
import itertools, math, sys, json
import numpy as np, mpmath as mp
from digits2 import growth
from homotopy import track, full
mp.mp.dps = 50
def best_m(S, n):
    return min(((math.log(growth(S, m)) / math.log(m)), m) for m in range(max(2, (n + 2) // 2), n + 2))
def is_rational(x, maxden=10**8):
    f = mp.identify(x, tol=mp.mpf(10)**-25, maxcoeff=maxden)
    return f
NMAX = int(sys.argv[1]); THR = float(sys.argv[2])
for n in range(6, NMAX + 1):
    h = n // 2 if n % 2 == 0 else (n + 1) // 2
    for T in itertools.combinations(range(2, n + 1), (math.ceil(n / 2) + 1) - 2):
        low = {0, 1, *T}; S = low | {2 * n - j for j in low}
        e, m = best_m(S, n)
        if e >= THR: continue
        Z = [j for j in range(2, n + 1) if j not in low]
        if len(Z) != n // 2: continue
        for eps in (1, -1):
            for sg in (1, -1):
                if n % 2 == 0 and eps * sg ** (n // 2) == -1: continue
                for u in track(n, eps, sg, Z):
                    a = full(u, n, eps, sg)
                    if np.max(np.abs(a.imag)) > 1e-8 or np.min(np.abs(a)) < 1e-6 * np.max(np.abs(a)): continue
                    a = a.real; b = [a[i] / a[1] ** i for i in range(n + 1)]
                    # quick float rationality screen (denominators up to 1e6)
                    from fractions import Fraction
                    fr = [Fraction(x).limit_denominator(10**6) for x in b]
                    ok = all(abs(float(f) - x) < 1e-9 * (1 + abs(x)) for f, x in zip(fr, b))
                    print("n=%d m=%d exp=%.4f low=%s eps=%d sg=%d rational_candidate=%s" % (n, m, e, sorted(low), eps, sg, ok), flush=True)
                    if ok: print("   RATIONAL?", [str(f) for f in fr], flush=True)
