"""All isolated solutions of the symmetric seed system via total-degree homotopy (quadratic equations).
Unknowns a_1..a_h (a_0=1), a_k for k>n/2 given by a_k = eps*sg^(n-k)*a_{n-k}. Equations: [x^j]p^2 = 0, j in Zlow."""
import numpy as np, json, sys, math
def full(u, n, eps, sg):
    h = n // 2; a = np.zeros(n + 1, complex); a[0] = 1; a[1:h+1] = u
    for k in range(h + 1, n + 1): a[k] = eps * sg ** (n - k) * a[n - k]
    return a
def F(u, n, eps, sg, Z):
    a = full(u, n, eps, sg); return np.array([np.dot(a[:j+1], a[j::-1]) for j in Z])
def J(u, n, eps, sg, Z):
    h = n // 2; a = full(u, n, eps, sg)
    # d a_k / d u_l
    D = np.zeros((n + 1, h), complex)
    for l in range(1, h + 1): D[l, l-1] = 1
    for k in range(h + 1, n + 1): D[k] = eps * sg ** (n - k) * D[n - k]
    return np.array([2 * sum(a[j - i] * D[i] for i in range(j + 1)) for j in Z])
def track(n, eps, sg, Z, seed=1):
    h = n // 2; assert len(Z) == h, (len(Z), h)
    rng = np.random.default_rng(seed); gamma = np.exp(2j * np.pi * rng.random())
    b = np.exp(2j * np.pi * rng.random(h))  # start system u_i^2 - b_i = 0
    G = lambda u: u ** 2 - b; GJ = lambda u: np.diag(2 * u)
    H = lambda u, t: (1 - t) * gamma * G(u) + t * F(u, n, eps, sg, Z)
    HJ = lambda u, t: (1 - t) * gamma * GJ(u) + t * J(u, n, eps, sg, Z)
    Ht = lambda u: F(u, n, eps, sg, Z) - gamma * G(u)
    sols = []
    for signs in range(2 ** h):
        u = np.array([math.sqrt(1) * np.sqrt(b[i]) * (1 if (signs >> i) & 1 else -1) for i in range(h)])
        t, dt, ok = 0.0, 0.01, True
        while t < 1:
            dt = min(dt, 1 - t)
            # Euler predictor
            try: du = np.linalg.solve(HJ(u, t), -Ht(u))
            except np.linalg.LinAlgError: ok = False; break
            u1 = u + dt * du; t1 = t + dt; conv = False
            for _ in range(6):
                try: c = np.linalg.solve(HJ(u1, t1), -H(u1, t1))
                except np.linalg.LinAlgError: break
                u1 = u1 + c
                if np.linalg.norm(c) < 1e-10 * (1 + np.linalg.norm(u1)): conv = True; break
            if conv: u, t = u1, t1; dt = min(dt * 1.5, 0.05)
            else:
                dt /= 2
                if dt < 1e-9: ok = False; break
            if np.linalg.norm(u) > 1e8: ok = False; break
        if ok:
            for _ in range(20):
                try: u = u + np.linalg.solve(J(u, n, eps, sg, Z), -F(u, n, eps, sg, Z))
                except np.linalg.LinAlgError: break
            if np.linalg.norm(F(u, n, eps, sg, Z)) < 1e-9: sols.append(u)
    uniq = []
    for s in sols:
        if not any(np.allclose(s, q, atol=1e-7) for q in uniq): uniq.append(s)
    return uniq
def run(n, low):
    Z = [j for j in range(2, n + 1) if j not in low]
    out = []
    for eps in (1, -1):
        for sg in (1, -1):
            if n % 2 == 0 and eps * sg ** (n // 2) == -1: continue
            if len(Z) != n // 2: continue
            sols = track(n, eps, sg, Z)
            nreal = 0
            for u in sols:
                a = full(u, n, eps, sg)
                real = np.max(np.abs(a.imag)) < 1e-8
                dense = np.min(np.abs(a)) > 1e-6 * np.max(np.abs(a))
                if real and dense:
                    nreal += 1; out.append((eps, sg, a.real))
            print("n=%d low=%s eps=%d sg=%d  #complex sols=%d  #real dense=%d" % (n, low, eps, sg, len(sols), nreal), flush=True)
    return out
if __name__ == "__main__":
    for n, low in json.loads(sys.argv[1]):
        for eps, sg, a in run(n, low):
            sq = np.convolve(a, a); nz = [i for i in range(2*n+1) if abs(sq[i]) > 1e-9 * np.abs(sq).max()]
            print("   REAL eps=%d sg=%d |p2|=%d supp=%s a=%s" % (eps, sg, len(nz), nz, json.dumps([round(float(v), 12) for v in a])))
