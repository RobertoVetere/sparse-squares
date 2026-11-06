"""Exhaustive enumeration of symmetric zero patterns with a square system (h-1 = ceil(n/2) unknowns/eqs).
S cap [0,n] = {0,1} u T, mirrored. Report patterns with lowest exponent log(lambda)/log(n+1)."""
import itertools, math, sys
from digits import growth
def run(n, top=8, slack=0):
    h = math.ceil(n / 2) + 1 + slack
    best = []
    for T in itertools.combinations(range(2, n + 1), h - 2):
        low = {0, 1, *T}
        S = low | {2 * n - j for j in low}
        lam = growth(S, n + 1)
        best.append((math.log(lam) / math.log(n + 1), lam, sorted(low)))
    best.sort()
    return best[:top]
if __name__ == "__main__":
    slack = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    for n in range(6, 25 if slack < 0 else 21):
        for e, lam, low in run(n, top=3, slack=slack):
            print(n, "exp=%.5f lam=%.4f" % (e, lam), "low=", low, flush=True)
