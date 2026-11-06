"""Search for solutions of OVERDETERMINED symmetric systems (one equation more than unknowns).
rho normalised to 1 via complex scaling; eps = +-1. Report any zero-residual solution."""
import numpy as np, math, itertools, sys
from scipy.optimize import least_squares
from digits import growth
def series(cvals, L, n):
    c = np.zeros(n + 1, complex); c[0] = 1; c[list(L)] = cvals
    a = np.zeros(n + 1, complex); a[0] = 1
    for j in range(1, n + 1):
        a[j] = (c[j] - np.dot(a[1:j], a[j-1:0:-1])) / 2
    return a
def res(x, L, n, eps):
    k = len(L); cv = x[:k] + 1j * x[k:]
    a = series(cv, L, n)
    m = (n + 1) // 2
    r = np.array([a[n - i] - eps * a[i] for i in range(m)])
    # scale-invariance guard: divide by norm of a
    r = r / (1 + np.abs(a).max())
    return np.concatenate([r.real, r.imag])
def search(n, low, tries=400, seed=0):
    L = sorted(set(low) - {0}); rng = np.random.default_rng(seed); found = []
    for eps in (1, -1):
        if n % 2 == 0 and eps == -1: continue
        for t in range(tries):
            x0 = rng.normal(size=2 * len(L)) * rng.choice([0.5, 2, 5])
            sol = least_squares(res, x0, args=(L, n, eps), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
            if np.linalg.norm(sol.fun) < 1e-11:
                cv = sol.x[:len(L)] + 1j * sol.x[len(L):]
                a = series(cv, L, n)
                if np.min(np.abs(a)) > 1e-6 * np.max(np.abs(a)):
                    found.append((eps, cv, a))
    return found
if __name__ == "__main__":
    import json
    targets = json.loads(sys.argv[1])
    for n, low in targets:
        S = set(low) | {2*n - j for j in low}
        print("n=%d low=%s exp_if_solvable=%.4f" % (n, low, math.log(growth(S, n+1))/math.log(n+1)), flush=True)
        f = search(n, low, tries=int(sys.argv[2]) if len(sys.argv) > 2 else 200)
        print("   solutions:", len(f), flush=True)
        for eps, cv, a in f[:5]:
            print("   eps", eps, "c", np.round(cv, 5), "a", np.round(a, 4))
