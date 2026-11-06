#!/usr/bin/env python3
"""Independent exact verifier for the Epoch "Dense Polynomials Whose Squares Are Sparse" problem.

Deliberately self-contained: pure Python integers, schoolbook squaring, no imports from
the search code. Usage:
    python3 verify.py candidate.json            # Epoch criterion (dense, Q(p^2) < sqrt(n))
    python3 verify.py --seed candidate.json     # report exponent log Q(p^2) / log T(p)
JSON: list of integers [a_0, ..., a_n] (constant term first).
"""
import json
import math
import sys
from fractions import Fraction


def square(coeffs):
    out = [0] * (2 * len(coeffs) - 1)
    for i, a in enumerate(coeffs):
        if a == 0:
            continue
        for j, b in enumerate(coeffs):
            if b:
                out[i + j] += a * b
    return out


def check_epoch(coeffs):
    assert all(type(c) is int for c in coeffs), "coefficients must be Python ints"
    n = len(coeffs) - 1
    dense = all(c != 0 for c in coeffs)
    q = sum(1 for c in square(coeffs) if c != 0)
    # strict inequality q < sqrt(n)  <=>  q*q < n  (q, n non-negative integers)
    ok = dense and n >= 1 and q * q < n
    return {"degree": n, "dense": dense, "terms_p": sum(1 for c in coeffs if c),
            "terms_p2": q, "sqrt_n": math.sqrt(n), "EPOCH_CRITERION_MET": ok}


def seed_report(coeffs):
    t = sum(1 for c in coeffs if c)
    q = sum(1 for c in square(coeffs) if c)
    return {"terms_p": t, "terms_p2": q, "dense": all(coeffs),
            "exponent_log_q_over_log_t": math.log(q) / math.log(t) if t > 1 else None}


if __name__ == "__main__":
    args = sys.argv[1:]
    seed = "--seed" in args
    args = [a for a in args if a != "--seed"]
    with open(args[0]) as fh:
        data = json.load(fh)
    if isinstance(data, dict):
        data = data["coefficients"]
    if any(isinstance(c, (str, Fraction)) for c in data):
        data = [int(c) for c in data]
    print(json.dumps(seed_report(data) if seed else check_epoch(data), indent=1))
