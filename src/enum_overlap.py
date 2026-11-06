import itertools, math, sys
from digits2 import growth
res = []
for n in range(4, 19):
    for slack in (0,):
        h = math.ceil(n/2) + 1 + slack
        for T in itertools.combinations(range(2, n+1), h - 2):
            low = {0, 1, *T}; S = low | {2*n - j for j in low}
            for m in range(max(2, (n+2)//2), n+2):
                g = growth(S, m)
                res.append((math.log(g)/math.log(m), n, m, g, sorted(low)))
res.sort()
seen=set()
for r in res[:40]:
    print("exp=%.5f n=%d m=%d lam=%.4f low=%s" % r)
