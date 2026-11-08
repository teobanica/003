import numpy as np
from itertools import combinations
from collections import Counter
from random import Random
from math import sqrt
def polar_check(U,r):
    A=U[:r,:r];B=U[:r,r:];C=U[r:,:r];D=U[r:,r:]
    left, singular, right = np.linalg.svd(A)
    V = left @ right
    R = (left * singular) @ left.conj().T
    left, singular, right = np.linalg.svd(D)
    W = left @ right
    T = (right.conj().T * singular) @ right
    K=np.linalg.inv(np.eye(r)+R)
    np.testing.assert_allclose(T,np.eye(len(D))-B.conj().T@K@B,atol=5e-13)
    np.testing.assert_allclose(W,D-C@V.conj().T@K@B,atol=5e-13)
def counts(E,I,J,p):
    out=Counter()
    for a,b in combinations(I,2):
        for c,d in combinations(J,2):
            val=(int(E[a,c])-int(E[b,c])-int(E[a,d])+int(E[b,d]))%p
            out[min(val,(-val)%p)]+=1
    return out
rng=Random(22)
for p in (7,11,13,17,19,23,29):
    m=(p-1)//2; Q=sorted({x*x%p for x in range(1,p)}); N=sorted(set(range(1,p))-set(Q))
    order=Q+[0]+N;E=np.outer(np.arange(p),np.arange(p))%p
    U=np.exp(2j*np.pi*E/p)/sqrt(p)
    U=U[np.ix_(order,order)]
    polar_check(U,m);polar_check(U,m+1)
    sN=sum(np.exp(2j*np.pi*np.array(N)/p))
    v=np.array([-sN]+[1]*m);v=v/np.linalg.norm(v);P=np.outer(v,v.conj())
    D=U[m:,m:];Dnew=D@(np.eye(m+1)+(1j-1)*P)
    Unew=U.copy();Unew[m:,m:]=Dnew
    np.testing.assert_allclose(Unew.conj().T@Unew,np.eye(p),atol=4e-13)
    sh=(Dnew-D)[1:,1:];np.testing.assert_allclose(sh,sh[0,0],atol=5e-13)
    assert np.max(abs(abs(Dnew)-1/sqrt(p)))>1e-5
    for I,J in [(Q,Q),(rng.sample(range(p),m),rng.sample(range(p),m))]:
        a=counts(E,I,J,p);d=counts(E,sorted(set(range(p))-set(I)),sorted(set(range(p))-set(J)),p)
        assert all(d[t]-a[t]==m*m for t in range(1,m+1))
        if p==19 and I==Q:assert set(a.values())=={144} and set(d.values())=={225}
print('Exact minor counts at seven primes; 144 and 225 at p=19. The accompanying polar checks use floating point.')
