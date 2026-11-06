import sys, json, math, numpy as np
from rational_strict import refine, ratcheck
from homotopy import track, full
from digits2 import growth
for n, low in json.loads(sys.argv[1]):
    Z = [j for j in range(2, n + 1) if j not in low]; S = set(low) | {2*n - j for j in low}
    e, m = min(((math.log(growth(S, m)) / math.log(m)), m) for m in range(max(2, (n + 2)//2), n + 2))
    tot = nr = 0
    for eps in (1, -1):
        for sg in (1, -1):
            if n % 2 == 0 and eps * sg ** (n // 2) == -1: continue
            if len(Z) != n // 2: continue
            for u in track(n, eps, sg, Z):
                a = full(u, n, eps, sg); tot += 1
                if not np.all(np.isfinite(a)) or np.min(np.abs(a)) < 1e-6 * np.max(np.abs(a)): continue
                A = refine(list(a), n, eps, sg, Z); b = [A[i] / A[1] ** i for i in range(n + 1)]
                rs = [ratcheck(x) for x in b]
                if all(r is not None for r in rs): nr += 1; print("RATIONAL", n, low, eps, sg, [str(r) for r in rs], flush=True)
    print("n=%d m=%d exp=%.4f  solutions checked=%d rational=%d" % (n, m, e, tot, nr), flush=True)
