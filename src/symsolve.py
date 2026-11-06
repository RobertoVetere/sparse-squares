"""Symmetric seed solver.
Unknowns: c_j (j in L = S cap [1,n]) with c_0 = 1. p = truncated sqrt series determined by the c's
(a_0..a_n). Symmetry ansatz a_{n-i} = eps * rho^i ... normalised rho = +-1:  a_{n-i} = eps * sgn^i * a_i.
Then p^2 vanishes on [0,n]\\S and, by symmetry, on the mirror. Solve by damped Newton from random complex starts,
keep real solutions with all a_i != 0."""
import numpy as np, mpmath as mp, random
def coeffs(c, n):
    # a = sqrt of C(x) = sum c_j x^j, a_0 = 1 (c_0 = 1), up to degree n
    a = [1.0 + 0j] + [0j] * n
    for j in range(1, n + 1):
        s = c.get(j, 0) - sum(a[i] * a[j - i] for i in range(1, j))
        a[j] = s / 2
    return a
def residual(x, L, n, eps, sg):
    c = dict(zip(L, x)); a = coeffs(c, n)
    return np.array([a[n - i] - eps * (sg ** i) * a[i] for i in range((n + 1) // 2)] if n % 2 else
                    [a[n - i] - eps * (sg ** i) * a[i] for i in range(n // 2)])
def newton(x, f, it=80):
    for _ in range(it):
        F = f(x); 
        if np.linalg.norm(F) < 1e-13: return x, True
        J = np.zeros((len(F), len(x)), complex); h = 1e-7
        for k in range(len(x)):
            e = np.zeros(len(x), complex); e[k] = h; J[:, k] = (f(x + e) - F) / h
        dx = np.linalg.lstsq(J, -F, rcond=None)[0]
        t = 1.0
        while t > 1e-4 and np.linalg.norm(f(x + t * dx)) > np.linalg.norm(F): t /= 2
        x = x + t * dx
    return x, np.linalg.norm(f(x)) < 1e-10
def solve(n, L, eps, sg, tries=300, scale=3.0, seed=0):
    rng = np.random.default_rng(seed); sols = []
    f = lambda x: residual(x, L, n, eps, sg)
    for _ in range(tries):
        x0 = (rng.normal(size=len(L)) + 1j * rng.normal(size=len(L))) * scale
        x, ok = newton(x0, f)
        if ok and np.all(np.abs(x.imag) < 1e-8 * (1 + np.abs(x.real))):
            xr = x.real
            a = np.array(coeffs(dict(zip(L, xr)), n)).real
            if np.min(np.abs(a)) > 1e-8 * np.max(np.abs(a)) and not any(np.allclose(xr, s, rtol=1e-6) for s in sols):
                sols.append(xr)
    return sols
