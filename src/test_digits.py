import random, math
from digits import *
for _ in range(300):
    b = random.randint(2, 7); S = set(random.sample(range(2*b-1), random.randint(1, 2*b-1)))
    k = random.randint(1, 4)
    assert count(S, b, k) == brute(S, b, k), (S, b, k)
print("automaton ok")
def mul(a,b_):
    r=[0]*(len(a)+len(b_)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b_): r[i+j]+=x*y
    return r
A=[1,2,-2,4,-10,50,125]
for al in [-110,-253]:
    B=[1,0,0,0,0,0,al]; p=mul(A,B); S={i for i,c in enumerate(mul(p,p)) if c}
    lam=growth(S,13); print(al, sorted(S), "lambda",lam, "exp", math.log(lam)/math.log(13))
