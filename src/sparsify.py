"""Structure-agnostic search: minimise a smooth surrogate of |supp(p^2)| over dense real p (a0=1, a1=1),
without prescribing the zero pattern. Continuation in delta; barrier keeps |a_i| away from 0 (relative).
Candidates are polished by Newton/LM on the detected zero set and accepted only if residual ~ 1e-13."""
import numpy as np, sys, time, json, math
from scipy.optimize import minimize, least_squares
def sq(a): return np.convolve(a, a)
def obj(x, n, delta, mu):
    a = np.concatenate([[1.0, 1.0], x]); c = sq(a); s = np.max(np.abs(a))
    cn = c / s ** 2
    f = np.sum(np.log(cn ** 2 + delta)) - mu * np.sum(np.log((a / s) ** 2 + 1e-300))
    # gradient
    g_c = 2 * cn / (cn ** 2 + delta) / s ** 2
    g_a = 2 * np.convolve(g_c, a[::-1])[n:2 * n + 1]  # d/da_k sum_j g_c[j] c_j = 2 sum_j g_c[j] a_{j-k}
    g_a = np.array([2 * np.dot(g_c[k:k + n + 1], a) for k in range(n + 1)])
    g_a -= mu * 2 / a
    return f, g_a[2:]
def polish(a, n, tol=1e-7):
    c = sq(a); Z = [j for j in range(2 * n + 1) if abs(c[j]) < tol * np.max(np.abs(c))]
    if not Z: return None
    res = lambda x: (lambda aa: sq(aa)[Z] / np.max(np.abs(aa)) ** 2)(np.concatenate([[1.0, 1.0], x]))
    r = least_squares(res, a[2:], xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
    b = np.concatenate([[1.0, 1.0], r.x])
    if np.linalg.norm(r.fun) > 1e-12 or np.min(np.abs(b)) < 1e-6 * np.max(np.abs(b)): return None
    c = sq(b); S = [j for j in range(2 * n + 1) if abs(c[j]) > 1e-9 * np.max(np.abs(c))]
    return b, S
def run(n, trials, seed, budget):
    rng = np.random.default_rng(seed); best = {}; t0 = time.time()
    for t in range(trials):
        if time.time() - t0 > budget: break
        x = rng.normal(size=n - 1) * rng.choice([0.5, 1, 3])
        for delta in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12]:
            r = minimize(obj, x, args=(n, delta, 1e-3), jac=True, method="L-BFGS-B", options=dict(maxiter=400))
            x = r.x
        a = np.concatenate([[1.0, 1.0], x]); p = polish(a, n)
        if p:
            b, S = p; k = len(S)
            if k not in best: best[k] = (b, S)
    return best
if __name__ == "__main__":
    n = int(sys.argv[1]); trials = int(sys.argv[2]); budget = float(sys.argv[3])
    best = run(n, trials, int(sys.argv[4]) if len(sys.argv) > 4 else 0, budget)
    for k in sorted(best)[:4]:
        b, S = best[k]; print("n=%d |p2|=%d (|p|=%d) S=%s" % (n, k, n + 1, S), flush=True)
