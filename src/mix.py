import itertools, math, numpy as np
from digits import transfer
pats = []
for n in range(4, 17):
    h = math.ceil(n/2) + 1
    for slack in (0, 1):
        for T in itertools.combinations(range(2, n+1), h - 2 + slack):
            low = {0, 1, *T}; S = low | {2*n - j for j in low}
            M, _ = transfer(S, n+1)
            rho = max(abs(np.linalg.eigvals(M)))
            pats.append((math.log(rho)/math.log(n+1), n, slack, tuple(sorted(low)), M))
pats.sort(key=lambda t: t[0])
print("patterns", len(pats)); 
# keep best 60 per n for mixing
byn = {}
for p in pats: byn.setdefault(p[1], []).append(p)
cand = [p for n in byn for p in byn[n][:25]]
best = []
for a, b in itertools.combinations_with_replacement(cand, 2):
    r = max(abs(np.linalg.eigvals(a[4] @ b[4])))
    e = math.log(r)/math.log((a[1]+1)*(b[1]+1))
    best.append((e, a[1], a[2], a[3], b[1], b[2], b[3], a[0], b[0]))
best.sort()
for x in best[:12]: print("exp=%.5f"%x[0], "n1=%d s=%d %s | n2=%d s=%d %s | single %.4f %.4f"%(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]))
