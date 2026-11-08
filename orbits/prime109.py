from itertools import product, combinations
from math import comb
p=109; g=36; n=6
square=[pow(g,j,p) for j in range(54)]; log={x:i for i,x in enumerate(square)}
assert len(log)==54
full_log={pow(n,j,p):j for j in range(108)}
def seed(t,r):
    row=[0]*54
    for x in range(1,p):
        y=(x+t)%p
        if y==0 or (full_log[y]-full_log[x])%4!=r:continue
        if r%2:
            row[log[x if x in log else y]]+=1
        elif x in log:
            row[log[y]]+=1;row[log[x]]-=1
    return row
def rem27(poly):
    a=poly[:]
    for k in range(len(a)-1,2,-1):
        c=a[k]
        a[k-1]=(a[k-1]-3*c)%27;a[k-2]=(a[k-2]-6*c)%27;a[k-3]=(a[k-3]+10*c)%27
    return tuple(z%27 for z in a[:3])
assert all(rem27(seed(t,r))==(0,0,0) for t,r in product([1,6],range(4)))
assert rem27([-1]+[0]*53+[1])==(0,0,0)
R=seed(6,1)
assert [sum(c*comb(j,i) for j,c in enumerate(R) if j>=i)%3 for i in range(4)]==[0,0,0,2]
assert sum((-1)**j*c for j,c in enumerate(R))%3==1
print('recurrence confirmed for all eight seed rows and full period')
for v,expected in [([0,0,1],[[2,0,2],[1,2,0],[2,0,0]]),([0,1,2],[[1,2,1],[0,2,1],[0,2,0]])]:
    # Direct differences of squares, with v the gap between the two position functions.
    # For A=(s+v)/2, D=(s-v)/2 the difference of squared gaps equals gap(s)*gap(v).
    M=[]
    for t,b in [(1,0),(1,1),(2,0)]:
        row=[]
        for a1,a2,k in [(1,0,0),(0,1,0),(0,0,1)]:
            S=lambda i:([0,a1,a2][i%3]+k*(i//3))%3
            A=lambda i:2*(S(i)+v[i%3])%3
            D=lambda i:2*(S(i)-v[i%3])%3
            pairs=[(i,j) for i,x in enumerate(square) for j,y in enumerate(square) if (y-x)%109==t and (j-i)%2==b]
            row.append(sum((A(j)-A(i))**2-(D(j)-D(i))**2 for i,j in pairs)%3)
        M.append(row)
    assert M==expected
    assert all(any(sum(x*y for x,y in zip(row,c))%3 for row in M) for c in product(range(3),repeat=3) if any(c))
    print('second moment rows',M,'force constants')
