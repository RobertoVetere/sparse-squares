import numpy as np, sys, json
from symsolve import solve, coeffs
targets = json.loads(sys.argv[1]); tries = int(sys.argv[2])
for n, low in targets:
    L = [j for j in low if j]
    for eps in (1, -1):
        for sg in (1, -1):
            if n % 2 == 0 and eps * sg ** (n // 2) == -1: continue
            for scale in (0.5, 2.0, 6.0):
                for x in solve(n, L, eps, sg, tries=tries, scale=scale, seed=int(scale*10)):
                    a = np.array(coeffs(dict(zip(L, x)), n)).real
                    sq = np.convolve(a, a); nz = [i for i in range(2*n+1) if abs(sq[i]) > 1e-9 * np.abs(sq).max()]
                    print(n, low, "eps", eps, "sg", sg, "|p2|", len(nz), "a=", json.dumps([float(v) for v in a]), flush=True)
