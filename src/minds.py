"""For tower base m with digits in [0,n] and K digits: for each N in [0,D], is the minimal-digit-sum
representation unique?  (Unique => coefficient nonzero after l-adic scaling, by valuation argument.)"""
import sys
def check(n, m, K):
    D = n * (m**K - 1) // (m - 1)
    # DP over positions: best[N] = (min digit sum, count of reps achieving it) for numbers representable with j digits
    best = {0: (0, 1)}
    for j in range(K):
        nb = {}
        w = m ** j
        for N, (ds, c) in best.items():
            for d in range(n + 1):
                M = N + d * w; v = (ds + d, c)
                if M not in nb or v[0] < nb[M][0]: nb[M] = v
                elif v[0] == nb[M][0]: nb[M] = (v[0], nb[M][1] + c)
        best = nb
    bad = [N for N in range(D + 1) if best[N][1] > 1]
    return D, len(bad), bad[:10]
for n, m in [(14, 12), (20, 17), (12, 11)]:
    for K in range(1, 5 if m < 15 else 4):
        print(n, m, K, check(n, m, K), flush=True)
