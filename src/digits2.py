"""General carry automaton: support of prod_i q(x^{m^i})^2 where supp(q^2)=S subset [0,D], any base m>=2.
Distinct sums sum s_i m^i. Deterministic subset construction over carries."""
import numpy as np
from itertools import product
def build(S, m):
    S = set(S); D = max(S); cmax = D // (m - 1) + 1
    start = frozenset([0]); states = {start: 0}; order = [start]; edges = []
    i = 0
    while i < len(order):
        Q = order[i]; i += 1
        for d in range(m):
            Qn = frozenset((s + c - d) // m for c in Q for s in S if (s + c - d) % m == 0 and s + c - d >= 0)
            if Qn:
                if Qn not in states: states[Qn] = len(order); order.append(Qn)
                edges.append((states[Q], states[Qn]))
    M = np.zeros((len(order), len(order)))
    for a, b in edges: M[a, b] += 1
    return M, order
def growth(S, m):
    M, _ = build(S, m)
    return max(abs(np.linalg.eigvals(M))) if len(M) else 0
def count(S, m, k):
    M, order = build(S, m); v = np.zeros(len(order)); v[0] = 1
    for _ in range(k): v = v @ M
    # final carry c contributes number c*m^k with its own digits; distinct (string, final carry value)
    tot = 0
    for j, Q in enumerate(order):
        # final carry values distinct numbers: count distinct values of carry
        tot += v[j] * len(Q)
    return int(round(tot))
def brute(S, m, k):
    return len({sum(s * m**i for i, s in enumerate(t)) for t in product(sorted(S), repeat=k)})
