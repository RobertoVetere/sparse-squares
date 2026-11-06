import itertools, math, sys, json
import numpy as np, mpmath as mp
from fractions import Fraction
from digits2 import growth
from homotopy import track, full
mp.mp.dps = 70
def refine(a, n, eps, sg, Z):
    h = n // 2; u = [mp.mpc(x) for x in a[1:h + 1]]
    def fullm(u):
        b = [mp.mpf(1)] + list(u)
        for k in range(h + 1, n + 1): b.append(eps * sg ** (n - k) * b[n - k])
        return b
    def F(u):
        b = fullm(u); return mp.matrix([mp.fsum(b[i] * b[j - i] for i in range(j + 1)) for j in Z])
    for _ in range(60):
        f0 = F(u); J = mp.matrix(len(Z), h); e = mp.mpf(10) ** -35
        for l in range(h):
            up = list(u); up[l] += e; J[:, l] = (F(up) - f0) / e
        du = mp.lu_solve(J, -f0); u = [u[i] + du[i] for i in range(h)]
        if mp.norm(du) < mp.mpf(10) ** -65: break
    return fullm(u)
def ratcheck(x):
    if abs(mp.im(x)) > mp.mpf(10)**-40: return None
    x = mp.re(x); fr = Fraction(str(mp.nstr(x, 65))).limit_denominator(10**12)
    return fr if abs(mp.mpf(fr.numerator) / fr.denominator - x) < mp.mpf(10) ** -45 else None
NMAX = int(sys.argv[1]); THR = float(sys.argv[2]); found = []
for n in range(6, NMAX + 1):
    for T in itertools.combinations(range(2, n + 1), (math.ceil(n / 2) + 1) - 2):
        low = {0, 1, *T}; S = low | {2 * n - j for j in low}
        e, m = min(((math.log(growth(S, m)) / math.log(m)), m) for m in range(max(2, (n + 2) // 2), n + 2))
        if e >= THR: continue
        Z = [j for j in range(2, n + 1) if j not in low]
        if len(Z) != n // 2: continue
        for eps in (1, -1):
            for sg in (1, -1):
                if n % 2 == 0 and eps * sg ** (n // 2) == -1: continue
                for u in track(n, eps, sg, Z):
                    a = full(u, n, eps, sg)
                    if not np.all(np.isfinite(a)) or np.min(np.abs(a)) < 1e-6 * np.max(np.abs(a)): continue
                    A = refine(list(a), n, eps, sg, Z)
                    # rescale x so a_1 = 1 (complex scaling allowed; rational => real seed up to x->tx)
                    b = [A[i] / A[1] ** i for i in range(n + 1)]
                    rs = [ratcheck(x) for x in b]
                    if all(r is not None for r in rs):
                        print("RATIONAL n=%d m=%d exp=%.5f low=%s eps=%d sg=%d b=%s" % (n, m, e, sorted(low), eps, sg, [str(r) for r in rs]), flush=True)
print("done")
