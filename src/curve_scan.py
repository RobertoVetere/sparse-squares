"""Rational-point scan on the slack-1 curve n=18, low=[0..6,15..18] (best m=15, exp 0.81137).
Normalisation a_0=1, a_1=1 (x-scaling). Symmetry: a_{n-i} = C rho^i a_i with C = rho^(-n/2).
For fixed rational rho the system in u=(a_2..a_9) is square & quadratic: solve all by total-degree homotopy."""
import numpy as np, math, sys, mpmath as mp
from fractions import Fraction
n = 18; low = set(range(7)) | {15, 16, 17, 18}; Z = [j for j in range(2, n + 1) if j not in low]; h = n // 2
def lin(rho):
    # a = A0 + B u
    A0 = np.zeros(n + 1, complex); B = np.zeros((n + 1, h - 1), complex)
    A0[0] = 1; A0[1] = 1
    for l in range(h - 1): B[l + 2, l] = 1
    C = rho ** (-n / 2)
    for k in range(h + 1, n + 1):
        i = n - k; A0[k] = C * rho ** i * A0[i]; B[k] = C * rho ** i * B[i]
    return A0, B
def solve_all(rho, seed=0):
    A0, B = lin(rho); d = h - 1; rng = np.random.default_rng(seed)
    def F(u):
        a = A0 + B @ u; return np.array([np.dot(a[:j + 1], a[j::-1]) for j in Z])
    def J(u):
        a = A0 + B @ u; return np.array([2 * sum(a[j - i] * B[i] for i in range(j + 1)) for j in Z])
    g = np.exp(2j * np.pi * rng.random()); b = np.exp(2j * np.pi * rng.random(d))
    sols = []
    for s in range(2 ** d):
        u = np.array([np.sqrt(b[i]) * (1 if (s >> i) & 1 else -1) for i in range(d)]); t = 0; dt = 0.02; ok = True
        while t < 1:
            dt = min(dt, 1 - t)
            HJ = lambda u, t: (1 - t) * g * np.diag(2 * u) + t * J(u)
            H = lambda u, t: (1 - t) * g * (u ** 2 - b) + t * F(u)
            try: du = np.linalg.solve(HJ(u, t), -(F(u) - g * (u ** 2 - b)))
            except np.linalg.LinAlgError: ok = False; break
            u1 = u + dt * du; t1 = t + dt; conv = False
            for _ in range(5):
                try: c = np.linalg.solve(HJ(u1, t1), -H(u1, t1))
                except np.linalg.LinAlgError: break
                u1 = u1 + c
                if np.linalg.norm(c) < 1e-10 * (1 + np.linalg.norm(u1)): conv = True; break
            if conv: u, t = u1, t1; dt = min(dt * 1.6, 0.08)
            else:
                dt /= 2
                if dt < 1e-8: ok = False; break
            if np.linalg.norm(u) > 1e7: ok = False; break
        if ok:
            for _ in range(15):
                try: u = u + np.linalg.solve(J(u), -F(u))
                except np.linalg.LinAlgError: break
            if np.linalg.norm(F(u)) < 1e-8: sols.append(A0 + B @ u)
    return sols
rhos = sorted({Fraction(p, q) for p in range(-9, 10) for q in range(1, 8) if p and math.gcd(p, q) == 1 and abs(Fraction(p,q)) not in (1,)}, key=lambda f: (abs(f.numerator) + f.denominator))
cnt = 0
for r in rhos[: int(sys.argv[1])]:
    if r < 0: continue  # C = rho^-9 real needs rho>0 for real... (negative also real: rho^-9 real) keep positive first
    sols = solve_all(float(r)); nreal = 0; hit = False
    for a in sols:
        if np.max(np.abs(a.imag)) < 1e-7 and np.min(np.abs(a)) > 1e-7 * np.max(np.abs(a)):
            nreal += 1
            fr = [Fraction(float(x.real)).limit_denominator(5000) for x in a]
            if all(abs(float(f) - x.real) < 1e-9 * (1 + abs(x.real)) for f, x in zip(fr, a)):
                hit = True; print("  RATIONAL CANDIDATE rho=%s a=%s" % (r, [str(f) for f in fr]), flush=True)
    print("rho=%s  #sols=%d  #real dense=%d" % (r, len(sols), nreal), flush=True)
