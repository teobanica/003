"""Exact prefix search by polynomial pseudo-remainders over the Eisenstein integers."""
from tangent import E, tangent
from functools import cache
from itertools import product
from collections import Counter
from math import gcd
import numpy as np
class W(tuple):
    def __new__(cls,a=0,b=0):return tuple.__new__(cls,(int(a),int(b)))
    def __add__(x,y):
        if type(y)!=W:y=W(y)
        return W(x[0]+y[0],x[1]+y[1])
    __radd__=__add__
    def __neg__(x):return W(-x[0],-x[1])
    def __sub__(x,y):return x+(-y)
    def __mul__(x,y):
        if type(y)!=W:y=W(y)
        a,b=x;c,d=y
        return W(a*c-b*d,a*d+b*c+b*d)
    __rmul__=__mul__
    def conj(x):return W(x[0]+x[1],-x[1])
def trim(v):
    v=list(v)
    while v and v[-1]==W():v.pop()
    if not v:return ()
    c=gcd(*(k for z in v for k in z));return tuple(W(z[0]//c,z[1]//c) for z in v)
def pgcd(a,b):
    while b:
        if len(b)==1:return (W(1),)
        r=a
        while len(r)>=len(b):
            d=len(r)-len(b);lead=r[-1]
            r=trim([b[-1]*v-(lead*b[j-d] if j>=d else W(0)) for j,v in enumerate(r)])
        a,b=b,r
    return a
six=[W(1),W(0,1),W(-1,1),W(-1),W(0,-1),W(1,-1)]
n=6; X,Y=E[:6,:6],E[:6,6:12]
k=[int(np.flatnonzero((x-y)%6==3)[0]) for x,y in zip(X,Y)]
def opts(i):return [(1,)*5] if i==0 else [(1,)+t for t in product([1,-1],repeat=4)]
@cache
def entries(i,s):
    signs=iter(opts(i)[s])
    return tuple((six[ix]-six[iy])*next(signs) if j!=k[i] else W() for j,(ix,iy) in enumerate(zip(X[i],Y[i])))
@cache
def coeff(i,si,j,sj):
    a,b=entries(i,si),entries(j,sj)
    return (2*six[X[i,k[i]]]*b[k[i]].conj(),2*six[X[j,k[j]]].conj()*a[k[j]],sum(x*y.conj() for x,y in zip(a,b)))
@cache
def norm(j,s):
    A,V,C=coeff(0,0,j,s)
    return trim((C*A.conj(),C*C.conj()+A*A.conj()-V*V.conj(),A*C.conj()))
@cache
def mix(i,si,j,sj):
    A,V,C=coeff(i,si,j,sj)
    Ai,Vi,Ci=[x.conj() for x in coeff(0,0,i,si)]
    Aj,Vj,Cj=coeff(0,0,j,sj)
    return trim((-A*Ai*Vj,C*Vi*Vj-A*Ci*Vj-V*Cj*Vi,-V*Aj*Vi))
levels=Counter();final=[]
def visit(sel,h):
    j=len(sel)
    if j==n:final.append((sel,h));return
    for v in range(16):
        g=pgcd(h,norm(j,v))
        for i in range(1,j):
            if len(g)==1:break
            g=pgcd(g,mix(i,sel[i],j,v))
        if len(g)==1:continue
        levels[j]+=1;visit(sel+[v],g)
visit([0],())
assert [levels[i] for i in range(1, 6)] == [16, 28, 4, 4, 4]
print('Surviving prefixes:', [levels[i] for i in range(1, 6)])
minors = {}
rows=[i for i in range(342) if i+1 not in {33,34,35,36,288,297,*range(331,343)}]
def detmod(A,q):
    a=A.copy()%q;v=1
    for k in range(len(a)):
        nz=np.flatnonzero(a[k:,k])
        if not len(nz):return 0
        r=k+nz[0]
        if r!=k:a[[r,k]]=a[[k,r]];v=-v
        pivot=int(a[k,k]);v=v*pivot%q
        fac=a[k+1:,k]*pow(pivot,-1,q)%q
        a[k+1:,k:]=(a[k+1:,k:]-fac[:,None]*a[k,k:])%q
    return v%q
for sel,P in final:
    assert len(P)==3 and P[1]==W() and P[0]==-P[2]
    print(sel,P)
    for e in (0,3):
        shifts=[e]
        for j in range(1,n):
            A,V,C=coeff(0,0,j,sel[j]); shifts.append(next(t for t in range(6) if A*six[e]+V*six[-t%6]+C==W()))
        H=E.copy()
        for i in range(n):
            signs=iter(opts(i)[sel[i]])
            for j in range(n):
                x,y=X[i,j],Y[i,j]
                if j==k[i]:x=(x+shifts[i])%6;y=(y+shifts[i])%6
                elif next(signs)<0:x,y=y,x
                H[i,j]=H[i+6,j+6]=x;H[i+6,j]=H[i,j+6]=y
        A=tangent(H)
        assert A.shape == (342, 324)
        det = detmod(A[rows], 101)
        assert det != 0
        if e == 0:
            minors[tuple(sel[1:])] = (tuple(shifts), det)
        print('first phase exponent', e, 'offsets', shifts, 'determinant mod 101', det)

# Printed representatives have first phase +1. Exponent 3 represents -1.
assert minors == {
    (0,0,0,0,0): ((0,0,0,0,0,0),19),
    (2,1,2,2,2): ((0,3,3,3,3,3),98),
    (9,1,5,9,14): ((0,0,3,0,0,3),46),
    (11,0,7,11,12): ((0,3,0,3,3,0),58),
}
