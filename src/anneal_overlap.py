import random, math, sys
from digits2 import growth
def score(low, n):
    S = set(low) | {2*n - j for j in low}
    best = (9, None)
    for m in range(max(2, (n + 2) // 2), n + 2):
        g = growth(S, m); e = math.log(g) / math.log(m)
        if e < best[0]: best = (e, m)
    return best
def anneal(n, slack, iters, seed):
    rng = random.Random(seed); h = math.ceil(n / 2) + 1 + slack
    pool = list(range(2, n + 1)); T = set(rng.sample(pool, h - 2))
    cur = score({0, 1} | T, n); best = (cur, sorted(T)); temp = 0.05
    for it in range(iters):
        o = rng.choice(sorted(T)); i = rng.choice([x for x in pool if x not in T])
        T2 = (T - {o}) | {i}; s = score({0, 1} | T2, n)
        if s[0] <= cur[0] or rng.random() < math.exp((cur[0] - s[0]) / temp):
            T, cur = T2, s
            if s[0] < best[0][0]: best = (s, sorted(T2))
        temp = max(0.002, temp * 0.997)
    return best
if __name__ == "__main__":
    for n in range(19, 33):
        for slack in (0,):
            r = min((anneal(n, slack, 1500, sd) for sd in range(3)), key=lambda t: t[0][0])
            print(n, "slack", slack, "exp=%.5f m=%d" % r[0], "low=", [0, 1] + r[1], flush=True)
