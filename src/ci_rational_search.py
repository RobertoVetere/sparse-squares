#!/usr/bin/env python3
"""Unattended search for RATIONAL seeds on positive-dimensional symmetric families (runs in GitHub Actions).

Family (n, low): p dense of degree n, a_0 = a_1 = 1 (x-scaling), symmetry
    a_{n-i} = sigma * sg^i * r^(2i-n) * a_i      (i < n/2),  r rational
zeros of p^2 at Z = [2..n] \\ low (mirrored zeros are automatic).
Unknowns: a_2..a_h and r (n even: rho=r enters as r^(i-n/2); n odd: r^(2i-n))  ->  slack = #unknowns - #equations >= 1 (curve/surface).
Sampling: fix r and (slack-1) of the a_i to random small-height rationals, solve the remaining square
quadratic system completely (total-degree homotopy), and test every solution for rationality at 60 digits.
A rational hit gives an integer seed; it is re-verified with exact integer arithmetic before being reported.
"""
import argparse, json, math, random, time, sys, os
from fractions import Fraction
import numpy as np
import mpmath as mp
from flint import fmpz_poly
sys.path.insert(0, os.path.dirname(__file__))
from digits2 import growth

mp.mp.dps = 60
THR = math.log(6) / math.log(9)  # integer record (Verdenius 1949)


def best_m(S, n):
    return min((math.log(growth(S, m)) / math.log(m), m) for m in range(max(2, (n + 2) // 2), n + 2))


def patterns(nmin, nmax, max_slack):
    out = []
    for n in range(nmin, nmax + 1):
        h = n // 2
        for r_ in range(1, n):
            for t in range(r_ + 2, n + 2):
                base = set(range(r_ + 1)) | set(range(t, n + 1))
                for extra in [None] + [x for x in range(r_ + 2, t - 1)]:
                    low = base | ({extra} if extra is not None else set())
                    Z = [j for j in range(2, n + 1) if j not in low]
                    slack = h - len(Z)  # unknowns a_2..a_h (h-1) + r (1)  minus equations
                    if not (1 <= slack <= max_slack):
                        continue
                    S = low | {2 * n - j for j in low}
                    e, m = best_m(S, n)
                    if e < THR:
                        out.append(dict(n=n, low=sorted(low), Z=Z, slack=slack, exp=e, m=m))
    out.sort(key=lambda d: d["exp"])
    return out


def lin(n, h, r, sigma, sg, fixed):
    """a = A0 + B u, u = free a_i (i in 2..h not fixed)."""
    free = [i for i in range(2, h + 1) if i not in fixed]
    A0 = np.zeros(n + 1, complex); B = np.zeros((n + 1, len(free)), complex)
    A0[0] = 1; A0[1] = 1
    for i, v in fixed.items(): A0[i] = float(v)
    for l, i in enumerate(free): B[i, l] = 1
    for k in range(h + 1, n + 1):
        i = n - k; f = sigma * sg ** i * float(r) ** ((i - n // 2) if n % 2 == 0 else (2 * i - n))
        A0[k] = f * A0[i]; B[k] = f * B[i]
    return A0, B, free


def homotopy(A0, B, Z, rng):
    d = B.shape[1]
    F = lambda u: np.array([np.dot((A0 + B @ u)[:j + 1], (A0 + B @ u)[j::-1]) for j in Z])
    def J(u):
        a = A0 + B @ u
        return np.array([2 * sum(a[j - i] * B[i] for i in range(j + 1)) for j in Z])
    g = np.exp(2j * np.pi * rng.random()); b = np.exp(2j * np.pi * rng.random(d))
    sols = []
    for s in range(2 ** d):
        u = np.array([np.sqrt(b[i]) * (1 if (s >> i) & 1 else -1) for i in range(d)]); t = 0.0; dt = 0.02; ok = True
        while t < 1:
            dt = min(dt, 1 - t)
            HJ = lambda u, t: (1 - t) * g * np.diag(2 * u) + t * J(u)
            H = lambda u, t: (1 - t) * g * (u ** 2 - b) + t * F(u)
            try:
                du = np.linalg.solve(HJ(u, t), -(F(u) - g * (u ** 2 - b)))
            except np.linalg.LinAlgError:
                ok = False; break
            u1 = u + dt * du; t1 = t + dt; conv = False
            for _ in range(5):
                try:
                    c = np.linalg.solve(HJ(u1, t1), -H(u1, t1))
                except np.linalg.LinAlgError:
                    break
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
            if np.linalg.norm(F(u)) < 1e-8: sols.append(u)
    return sols


def rational_check(u0, n, h, r, sigma, sg, fixed, free, Z):
    """Refine at 60 digits; return exact Fractions for all a_i if they are rational, else None."""
    u = [mp.mpc(x) for x in u0]
    def full(u):
        a = [mp.mpf(1), mp.mpf(1)] + [None] * (n - 1)
        for i, v in fixed.items(): a[i] = mp.mpf(v.numerator) / v.denominator
        for l, i in enumerate(free): a[i] = u[l]
        rr = mp.mpf(r.numerator) / r.denominator
        for k in range(h + 1, n + 1):
            i = n - k; a[k] = sigma * sg ** i * rr ** ((i - n // 2) if n % 2 == 0 else (2 * i - n)) * a[i]
        return a
    def F(u):
        a = full(u); return mp.matrix([mp.fsum(a[i] * a[j - i] for i in range(j + 1)) for j in Z])
    for _ in range(40):
        f0 = F(u); Jm = mp.matrix(len(Z), len(free)); e = mp.mpf(10) ** -30
        for l in range(len(free)):
            up = list(u); up[l] += e; Jm[:, l] = (F(up) - f0) / e
        try: du = mp.lu_solve(Jm, -f0)
        except ZeroDivisionError: return None
        u = [u[i] + du[i] for i in range(len(free))]
        if mp.norm(du) < mp.mpf(10) ** -55: break
    out = []
    for x in full(u):
        if abs(mp.im(x)) > mp.mpf(10) ** -40: return None
        fr = Fraction(mp.nstr(mp.re(x), 58)).limit_denominator(10 ** 14)
        if abs(mp.mpf(fr.numerator) / fr.denominator - mp.re(x)) > mp.mpf(10) ** -42: return None
        out.append(fr)
    return out


def exact_verify(fr, low, n):
    den = 1
    for f in fr: den = den * f.denominator // math.gcd(den, f.denominator)
    p = [int(f * den) for f in fr]
    g = 0
    for c in p: g = math.gcd(g, c)
    p = [c // g for c in p]
    sq = [int(c) for c in (fmpz_poly(p) ** 2).coeffs()]
    S = [j for j, c in enumerate(sq) if c]
    expected = sorted(set(low) | {2 * n - j for j in low})
    return p, S, all(p) and set(S) <= set(expected)


def small_rational(rng, H):
    while True:
        a, b = rng.randint(-H, H), rng.randint(1, H)
        if a and math.gcd(a, b) == 1: return Fraction(a, b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--budget", type=float, default=600); ap.add_argument("--nmax", type=int, default=24)
    ap.add_argument("--maxpaths", type=int, default=2048); ap.add_argument("--out", default="out.jsonl")
    args = ap.parse_args()
    rng = random.Random(1000 + args.shard); nrng = np.random.default_rng(args.shard)
    pats = [p for p in patterns(10, args.nmax, 3) if 2 ** (p["n"] // 2 - 1 - (p["slack"] - 1)) <= args.maxpaths]
    mine = pats[args.shard::args.nshards] or pats
    t0 = time.time(); stats = dict(shard=args.shard, patterns=len(mine), samples=0, solutions=0, real_dense=0, rational_hits=0)
    with open(args.out, "w") as fh:
        fh.write(json.dumps({"type": "patterns", "shard": args.shard, "list": [(p["n"], p["low"], p["slack"], round(p["exp"], 5), p["m"]) for p in mine]}) + "\n")
        k = 0
        while time.time() - t0 < args.budget:
            P = mine[k % len(mine)]; k += 1
            n, h, Z = P["n"], P["n"] // 2, P["Z"]
            H = rng.choice([3, 5, 8, 12])
            r = small_rational(rng, H); sg = rng.choice([1, -1])
            sigma = 1 if n % 2 == 0 else rng.choice([1, -1])
            idx = rng.sample(range(2, h + 1), P["slack"] - 1)
            fixed = {i: small_rational(rng, H) for i in idx}
            A0, B, free = lin(n, h, r, sigma, sg, fixed)
            if len(free) != len(Z): continue
            sols = homotopy(A0, B, Z, nrng); stats["samples"] += 1; stats["solutions"] += len(sols)
            for u in sols:
                a = A0 + B @ u
                if np.max(np.abs(a.imag)) < 1e-7 and np.min(np.abs(a)) > 1e-7 * np.max(np.abs(a)):
                    stats["real_dense"] += 1
                    fr = rational_check(list(u), n, h, r, sigma, sg, fixed, free, Z)
                    if fr is not None:
                        p, S, ok = exact_verify(fr, P["low"], n)
                        e, m = best_m(S, n)
                        rec = {"type": "RATIONAL_HIT", "n": n, "low": P["low"], "r": str(r), "sg": sg, "sigma": sigma,
                               "fixed": {str(i): str(v) for i, v in fixed.items()}, "p": p, "supp_p2": S,
                               "exact_ok": ok, "exp": e, "m": m, "beats_integer_record": e < THR}
                        stats["rational_hits"] += 1; fh.write(json.dumps(rec) + "\n"); fh.flush()
                        print("RATIONAL HIT", rec, flush=True)
            if stats["samples"] % 5 == 0:
                print(json.dumps(stats), "elapsed=%.0f" % (time.time() - t0), flush=True)
        stats["elapsed"] = time.time() - t0
        fh.write(json.dumps({"type": "stats", **stats}) + "\n")
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
