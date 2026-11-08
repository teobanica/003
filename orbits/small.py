import subprocess
from collections import defaultdict
import sys
from pathlib import Path
from sympy import primitive_root
from itertools import product
# The runner passes the two compiled enumerators as absolute paths.
pair_sums, differences = map(Path, sys.argv[1:])
expected = {5:(1,1,1,0,0),7:(2,2,2,0,0),11:(4,4,4,0,0),13:(12,6,2,0,0),
            17:(56,28,4,0,0),19:(66,66,6,54,6),23:(230,230,10,110,10),
            29:(396,198,6,672,6),31:(552,552,8,120,8)}
for p in (5,7,11,13,17,19,23,29,31):
    m=(p-1)//2
    def enumerate_words(program):
        lines = subprocess.check_output([program, str(p)]).splitlines()
        return sorted(tuple(map(int, line.split())) for line in lines)
    rows = enumerate_words(pair_sums)
    assert rows == enumerate_words(differences)
    assert len(set(rows)) == len(rows)
    n=int(primitive_root(p))
    def sig(word):
        return tuple(frozenset((word[(j+s)%m]-x)%p for j,x in enumerate(word)) for s in range(1,m))
    def geo(u):
        r=u[1]*pow(u[0],-1,p)%p
        return all(u[(i+1)%m]==x*r%p for i,x in enumerate(u))
    sq=defaultdict(list); ns=defaultdict(list)
    for a in rows:
        sq[sig(a)].append(a)
        b=tuple(x*n%p for x in a);ns[sig(b)].append(b)
    count=tests=accepted=0
    for k,A in sq.items():
        B=ns[tuple(frozenset(range(1,p))-d for d in k)]
        if B: assert all(map(geo,A+B))
        if not B:continue
        count+=len(A)
        if p<19:continue
        for a,d,b,c in product(A,A,B,B):
            for j in range(m):
                t=d[j:]+d[:j];tests+=1
                if all(len({(a[(i+s)%m]-c[i])%p for i in range(m)}|{(b[(i+s)%m]-t[i])%p for i in range(m)})==2*m for s in range(m)):
                    assert all(b[i]==c[i]==n*x%p and t[i]==n*n*x%p for i,x in enumerate(a));accepted+=1
    result = (len(rows),len(sq),count,tests,accepted)
    assert result == expected[p], (p, result)
    print(p, *result, flush=True)
