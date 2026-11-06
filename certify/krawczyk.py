#!/usr/bin/env python3
"""Rigorous certification (interval arithmetic, mpmath.iv) of a real symmetric seed.

System: p(x) = sum a_i x^i, a_0 = 1, unknowns u = (a_1..a_h), h = n//2,
        a_k = eps * sg^(n-k) * a_{n-k} for k > h  (requires eps*sg^(n/2) = 1 when n even),
        equations: [x^j] p(x)^2 = 0 for j in Z = [2..n] minus `low`.
Krawczyk: if K(X) is contained in the interior of box X, the system has exactly one solution in X.
Then the exact solution's coefficients lie in the printed intervals; we check every a_i excludes 0
and that [x^j]p^2 excludes 0 for j in S (so supp(p^2) = S exactly; zeros on Z hold exactly by construction).
Usage: krawczyk.py seed.json   (fields n, low, eps, sg, a)
"""
import json, sys
import mpmath as mp
from mpmath import iv
mp.mp.dps = 110; iv.dps = 110

def full(u, n, eps, sg, one):
    h = n // 2; a = [one] + list(u)
    for k in range(h + 1, n + 1): a.append(eps * sg ** (n - k) * a[n - k])
    return a

def main(path):
    d = json.load(open(path)); n, low, eps, sg = d["n"], set(d["low"]), d["eps"], d["sg"]
    h = n // 2; Z = [j for j in range(2, n + 1) if j not in low]
    assert len(Z) == h, "not a square system"
    if n % 2 == 0: assert eps * sg ** (n // 2) == 1
    # dependency matrix D[k][l] = d a_k / d u_l  (constant)
    D = [[0] * h for _ in range(n + 1)]
    for l in range(h): D[l + 1][l] = 1
    for k in range(h + 1, n + 1): D[k] = [eps * sg ** (n - k) * x for x in D[n - k]]
    def F(u, one):
        a = full(u, n, eps, sg, one); return [sum((a[i] * a[j - i] for i in range(j + 1)), 0 * one) for j in Z]
    def J(u, one):
        a = full(u, n, eps, sg, one)
        return [[sum((2 * a[j - i] * D[i][l] for i in range(j + 1) if D[i][l]), 0 * one) for l in range(h)] for j in Z]
    # high-precision Newton refinement of the midpoint
    u = [mp.mpf(x) for x in d["a"][1:h + 1]]
    for _ in range(100):
        Fm = mp.matrix(F(u, mp.mpf(1))); Jm = mp.matrix(J(u, mp.mpf(1)))
        du = mp.lu_solve(Jm, -Fm); u = [u[i] + du[i] for i in range(h)]
        if mp.norm(du) < mp.mpf(10) ** -100: break
    Y = mp.inverse(mp.matrix(J(u, mp.mpf(1))))
    delta = mp.mpf(10) ** -60
    X = [iv.mpf([x - delta, x + delta]) for x in u]
    one = iv.mpf(1)
    Fx = F([iv.mpf(x) for x in u], one); JX = J(X, one)
    Kb = []
    for i in range(h):
        s = iv.mpf(u[i]) - sum((iv.mpf(Y[i, k]) * Fx[k] for k in range(h)), iv.mpf(0))
        for l in range(h):
            c = (1 if i == l else 0) - sum((iv.mpf(Y[i, k]) * JX[k][l] for k in range(h)), iv.mpf(0))
            s += c * (X[l] - iv.mpf(u[l]))
        Kb.append(s)
    inside = all((Kb[i].a > X[i].a) and (Kb[i].b < X[i].b) for i in range(h))
    a = full(X, n, eps, sg, one)
    nonzero_a = all(not (x.a <= 0 <= x.b) for x in a)
    sq = [sum((a[i] * a[j - i] for i in range(max(0, j - n), min(j, n) + 1)), iv.mpf(0)) for j in range(2 * n + 1)]
    S = sorted(set(low) | {2 * n - j for j in low})
    support_ok = all(not (sq[j].a <= 0 <= sq[j].b) for j in S)
    # zeros off S: positions in Z hold exactly; mirrored positions hold by symmetry (checked as containment of 0)
    zeros_contain0 = all(sq[j].a <= 0 <= sq[j].b for j in range(2 * n + 1) if j not in S)
    out = {"n": n, "krawczyk_inclusion": bool(inside), "all_coefficients_nonzero": bool(nonzero_a),
           "square_support_exactly_S": bool(support_ok and zeros_contain0), "S": S, "|S|": len(S),
           "box_radius": "1e-60",
           "a_mid_30digits": [mp.nstr(mp.mpf(x.mid), 30) for x in a]}
    print(json.dumps(out, indent=1))
    return out

if __name__ == "__main__":
    main(sys.argv[1])
