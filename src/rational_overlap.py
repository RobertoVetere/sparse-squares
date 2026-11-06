"""Known RATIONAL dense seeds under overlapping bases (single seed and 2-seed mixed radix)."""
import math, numpy as np
from flint import fmpz_poly
from digits2 import growth
def supp(p):
    q=fmpz_poly(p); return [i for i,c in enumerate((q*q).coeffs()) if c!=0]
def mul(a,b):
    return (fmpz_poly(a)*fmpz_poly(b)).coeffs()
A=[1,2,-2,4,-10,50,125]
seeds={"V8":[1,-2,-2,-4,-10,-28,-98,-686,2401],"A6":A}
for al in [-110,-253]: seeds["CD12_%d"%al]=[int(c) for c in mul(A,[1,0,0,0,0,0,al])]
seeds["CD12_15625"]=[int(c) for c in mul(A,[1,0,0,0,0,0,15625])]
for k,p in seeds.items():
    n=len(p)-1; S=supp(p)
    best=min((math.log(growth(S,m))/math.log(m),m) for m in range(max(2,(n+2)//2),n+2))
    print(k,"n=%d |p2|=%d"%(n,len(S)),"best exp=%.5f at m=%d"%best)
