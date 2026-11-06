"""For base b (deg p = b-1), anneal over zero sets Z of p^2 (|Z| <= m) minimizing growth lambda.
Positions 0,1,2b-3,2b-2 are always nonzero for dense p."""
import random, math, sys
from digits import growth
def best(b, m, iters=20000, seed=0):
    rng = random.Random(seed)
    free = list(range(2, 2*b-3))
    m = min(m, len(free))
    Z = set(rng.sample(free, m))
    f = lambda Z: growth(set(range(2*b-1)) - Z, b)
    cur = f(Z); bestv, bestZ = cur, set(Z); T = 0.5
    for it in range(iters):
        if len(Z) >= len(free): break
        out = rng.choice(sorted(Z)); inn = rng.choice([x for x in free if x not in Z])
        Z2 = (Z - {out}) | {inn}; v = f(Z2)
        if v <= cur or rng.random() < math.exp((cur - v) / T):
            Z, cur = Z2, v
            if v < bestv: bestv, bestZ = v, set(Z2)
        T = max(0.01, T * 0.9995)
    return bestv, sorted(bestZ)
if __name__ == "__main__":
    for b in range(5, 22):
        for extra in (-1, 0, 1, 2):
            m = b - 2 + extra   # n-1 = b-2 equations generic count
            v, Z = best(b, m, iters=4000)
            print(b, "|Z|=%d" % m, "lam=%.4f exp=%.4f" % (v, math.log(v)/math.log(b)), "S=", sorted(set(range(2*b-1))-set(Z)), flush=True)
