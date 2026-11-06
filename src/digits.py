"""Growth rate of support of P_k^2 where P_k = prod_i p(x^{b^i}), b = deg p + 1.
Support of P_k^2 subset of { sum s_i b^i : s_i in S }, S = supp(p^2) subset [0, 2b-2].
Exact count via determinized carry automaton (states = sets of carries in {0,1})."""
import numpy as np
from itertools import product

def transfer(S, b):
    states = [frozenset([0]), frozenset([1]), frozenset([0, 1])]
    idx = {s: i for i, s in enumerate(states)}
    M = np.zeros((3, 3))
    for Q in states:
        for d in range(b):
            Qn = frozenset(cn for c in Q for cn in (0, 1) if (d + b * cn - c) in S)
            if Qn:
                M[idx[Q], idx[Qn]] += 1
    return M, idx

def count(S, b, k):
    S = set(S); M, idx = transfer(S, b)
    v = np.zeros(3); v[idx[frozenset([0])]] = 1
    for _ in range(k): v = v @ M
    return int(round(v[0] * 1 + v[1] * 1 + v[2] * 2))

def growth(S, b):
    M, _ = transfer(set(S), b)
    return max(abs(np.linalg.eigvals(M)))

def brute(S, b, k):
    return len({sum(s * b**i for i, s in enumerate(t)) for t in product(sorted(S), repeat=k)})
