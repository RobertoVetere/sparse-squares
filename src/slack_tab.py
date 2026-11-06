import itertools, math
from digits2 import growth
res = []
for n in range(6, 19):
    h = math.ceil(n/2) + 2
    for T in itertools.combinations(range(2, n+1), h - 2):
        low = {0, 1, *T}; S = low | {2*n - j for j in low}
        e, m = min(((math.log(growth(S, m))/math.log(m)), m) for m in range(max(2,(n+2)//2), n+2))
        res.append((e, n, m, sorted(low)))
res.sort()
for r in res[:15]: print("exp=%.5f n=%d m=%d low=%s" % r)
