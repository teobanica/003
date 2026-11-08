from itertools import permutations
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from sympy import factorint
from sympy.polys.matrices import DomainMatrix
from sympy import GF
def census(p,n):
    m=(p-1)//2;h=m//2
    F=[pow(n,4*i,p) for i in range(h)];good=[]; sig=[]
    for tail in permutations(F[1:]):
        u=(1,)+tail
        for j in range(1,h//2+1):
            if len({(u[(a+j)%h]-x)%p for a,x in enumerate(u)})!=h:break
        else:good.append(u)
    seq=[np.array(good)*pow(n,j,p)%p for j in range(4)]
    for X in seq:
        vals=np.stack([(np.roll(X,-s,axis=1)-X)%p for s in range(1,h//2+1)],axis=1)
        sig.append(np.sum(1<<vals,axis=2))
    alls=(1<<p)-2; halved=defaultdict(list);npairs=0;nquads=0;sol=[]
    for k in [0,1]:
        for i,left in enumerate(sig[k]):
            occ=np.flatnonzero(~np.any(left&sig[k+2],axis=1))
            for j in occ:
                key=tuple((left|sig[k+2][j])^(alls if k else 0))
                if k==0:halved[key].append((i,j));npairs+=1;continue
                for A,C in halved.get(key,[]):
                    nquads+=1
                    v=[0]*m;y=[0]*m;v[::2]=seq[0][A];y[::2]=seq[1][i]
                    for s in range(h):
                        v[1::2]=np.roll(seq[2][C],s)
                        for t in range(h):
                            y[1::2]=np.roll(seq[3][j],t)
                            if all(len({(u[(k+z)%m]-u[k])%p for u in [v,y] for k in range(m)})==p-1 for z in range(1,m)):
                                sol.append((v.copy(),y.copy()))
    assert all(all(len({u[(i+1)%m]*pow(int(u[i]),-1,p)%p for i in range(m)})==1 for u in pair) for pair in sol)
    return len(good),npairs,len(halved),nquads,len(sol)
def ranks(p,n):
    m=(p-1)//2;h=m//2; idx={pow(n,i,p):i for i in range(2*m)}
    Q=[x for x in range(1,p) if idx[x]%2==0]; N=set(idx)-set(Q);q=set(Q)
    # Build each matrix one row at a time directly from its defining equations.
    B=[]
    for t in range(1,m+1):
        row=[]
        for z in range(1,m+1):
            S=q if z in q else N
            row.append(sum(int((u-t)%p in S) for u in (z,-z))%2)
        B.append(row)
    for S in (q,N):B.append([int(i in S) for i in range(1,m+1)])
    rank=lambda L,l: DomainMatrix.from_list(L,GF(l)).to_sparse().rank()
    r=rank(B,2);L=[]
    for t in range(1,p):
        rows=[dict.fromkeys(Q,0) for _ in range(4)]
        for x in range(1,p):
            y=(x+t)%p
            if y==0:continue
            j=(idx[y]-idx[x])%4;z=rows[j]
            if j%2:z[x if x in q else y]+=1
            elif x in q:z[x]-=1;z[y]+=1
        L.extend([list(z.values()) for z in rows])
    out={l:rank(L,l) for l in factorint(h)}
    return r,out

if __name__ == '__main__':
    expected_ranks = {
        (37,2):(18,{3:17}), (41,6):(20,{2:18,5:19}), (53,2):(26,{13:25}),
        (61,2):(30,{3:29,5:29}), (89,3):(44,{2:43,11:43}), (149,2):(74,{37:73}),
        (173,2):(86,{43:85}), (197,2):(98,{7:97}), (233,3):(116,{2:115,29:115}),
    }
    for (p, n), expected in expected_ranks.items():
        result = ranks(p, n)
        assert result == expected, (p, result)
        print(p, 'ranks', result, flush=True)
    for p, n, expected in [(37, 2, (9378,888,444,12,12)), (41, 6, (4620,520,130,32,16))]:
        result = census(p, n)
        assert result == expected, (p, result)
        print(p, 'half-word counts', result, flush=True)
