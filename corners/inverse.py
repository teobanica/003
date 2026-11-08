from itertools import permutations
from pathlib import Path
from collections import Counter
from time import monotonic
def check(p):
    m=(p-1)//2
    Q=tuple(sorted({x*x%p for x in range(1,p)})); A=set(Q)
    N=tuple(v for v in range(1,p) if v not in A)
    chi=lambda d: int(d%p in A)
    b=Q[-1];dom=[a for a in Q if a!=b]
    # Enumerate every injective character-preserving map that fixes 1 on a one-vertex deletion.
    extend=[]
    def visit(P,T):
        if len(P)==len(dom):
            t,=A-set(P)
            mapping=dict(zip(dom,P));mapping[b]=t
            assert all(chi(u-v)==chi(mapping[u]-mapping[v]) for u in Q for v in Q if u!=v)
            extend.append(tuple(mapping[x] for x in Q));return
        a=dom[len(P)]
        for c in T:
            if a==1 and c!=1:continue
            if all(chi(a-u)==chi(c-v) for u,v in zip(dom,P)):visit(P+[c],T-{c})
    visit([],A);expected={tuple(Q)}
    if p%4==1:expected.add(tuple(pow(a,-1,p) for a in Q))
    assert set(extend)==expected
    if p%4==1 and m>=6:
        for y in N:
            for z in N:
                if y!=z:
                    entries=[(y*a-z*pow(a,-1,p))%p for a in dom]
                    assert len(entries)!=len(set(entries))
    print(p,'vertex extensions',len(extend),flush=True)
    if p not in [13,17]:return
    # Independently enumerate coset-exchanging permutations orthogonal to all retained rows.
    # In characteristic p these differences must be permutations, hence collisions prune immediately.
    xs=[u for pair in zip(Q,N) for u in pair]
    f={x:tuple((a*x if x in A else pow(a,-1,p)*x)%p for a in dom) for x in xs}
    used=[{0} for _ in dom]
    count=0
    def search(k,freeQ,freeN):
        nonlocal count
        if k==len(xs):count+=1;return
        x=xs[k]; target=freeN if x in A else freeQ
        for y in sorted(target):
            v=[(y-h)%p for h in f[x]]
            if any(h in t for h,t in zip(v,used)):continue
            for h,t in zip(v,used):t.add(h)
            target.remove(y);search(k+1,freeQ,freeN);target.add(y)
            for h,t in zip(v,used):t.remove(h)
    search(0,set(Q),set(N))
    assert count==0
    print(p,'coset-exchanging rows in the deleted inverse problem',count,flush=True)
for p in [5,7,11,13,17,19,23,29,31,37,41,43]:check(p)
