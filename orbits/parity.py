import numpy as np
from itertools import product
from even import ranks
from sympy import primitive_root

def reduce(A, n):
    A=A.copy(); piv=[]; r=0
    for j in range(n):
        ix=np.flatnonzero(A[r:,j])
        if not ix.size: continue
        k=r+ix[0]
        A[[r,k]]=A[[k,r]]
        take=A[:,j].copy();take[r]=False
        A[take]^=A[r];piv.append(j);r+=1
    return A,piv

def check(p):
    m=(p-1)//2;h=m//2;y=np.arange(1,m+1)
    Q=np.array([i in {j*j%p for j in range(1,p)} for i in range(p)])
    # Literal formula for M; zero is in neither class.
    def same(x,u):
        v=(x%p)
        return (v!=0)&(Q[v]==Q[u])
    M=np.stack([same(y-t,y)^same(-y-t,y) for t in y]+[Q[y],~Q[y]])
    aug=np.column_stack((M,np.full(m+2,bool(h%2))));R,pi=reduce(aug,m)
    free=sorted(set(range(m))-set(pi)); assert not R[len(pi):,-1].any()
    N=(~Q).copy();N[0]=False
    fields=np.arange(1,p)
    masks=[same(fields+t,fields) for t in y]
    xidx=[fields[mask] for mask in masks];yidx=[(i+t)%p for i,t in zip(xidx,y)]
    reps=np.minimum(np.arange(p),p-np.arange(p))
    neg=np.arange(p)>m
    flips=(np.arange(m+1)[:,None]==reps[None,:]);flips[0,:]=False
    survivors=[]
    for bits in product((False,True),repeat=len(free)):
        z=np.zeros(m,bool);z[free]=bits
        z[pi]=R[:len(pi),-1]^(R[:len(pi),free]@np.array(bits,np.int64)%2).astype(bool)
        sig=1-2*(flips^(np.r_[False,z][reps]&neg)).astype(np.int16)
        F=np.stack([1+(sig[:,xx]*sig[:,yy]).sum(axis=1) for xx,yy in zip(xidx,yidx)]+[sig[:,Q].sum(axis=1),sig[:,N].sum(axis=1)])
        base=F[:,0]
        if (base%4).any():continue
        diff=F[:,1:]-base[:,None];assert not (diff%4).any()
        A=np.column_stack(((diff//4)%2,(base//4)%2)).astype(bool)
        RA,piv=reduce(A,m)
        if not RA[len(piv):,-1].any():survivors.append(z)
    assert len(survivors)==1 and np.all(survivors[0]==h%2)
    matrix_rank, prime_ranks = ranks(p, int(primitive_root(p)))
    assert matrix_rank == len(pi)
    observed = (len(pi), 2**len(free), len(survivors), prime_ranks)
    expected = {73: (32,16,1,{2:35,3:35}), 97: (40,256,1,{2:47,3:47}),
                101: (46,16,1,{5:49}), 109: (52,4,1,{3:51})}
    assert observed == expected[p], (p, observed)
    print(p, observed, flush=True)
for p in [73,97,101,109]:check(p)
