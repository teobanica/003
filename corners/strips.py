from sympy import QQ
from time import monotonic

def solve(M, R, K):
    n = len(M)
    a = [row[:] + list(rhs) for row, rhs in zip(M, R)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j] != K.zero)
        a[pivot], a[j] = a[j], a[pivot]
        d = a[j][j]
        a[j] = [v/d for v in a[j]]
        for i in range(n):
            if i != j and a[i][j] != K.zero:
                d = a[i][j]
                a[i] = [x-d*y for x,y in zip(a[i],a[j])]
    return [row[n:] for row in a]

def run(p, rows, cols):
    start = monotonic()
    K = QQ.cyclotomic_field(p)
    zero, one, w = K.zero, K.one, K.unit
    roots = [w**j for j in range(p)]
    def gal(a, n):
        h = zero
        for c in a.to_list():
            h = h*roots[n % p] + K.convert(c)
        return h
    add = lambda values: sum(values, zero)
    R = [k for k in range(p) if k not in rows][::-1]
    D = [k for k in range(p) if k not in cols][::-1]
    # Independent permutations and phases outside the known Fourier intersection.
    def H(i,j):
        return roots[(i*j + (i*i+3 if i in R else 0) + (2*j+1 if j in D else 0))%p]
    P = [[H(i,j) for j in cols] for i in rows]
    B = [[H(i,j) for j in D] for i in rows]
    C = [[H(i,j) for j in cols] for i in R]
    r,s,n,l = len(rows),len(cols),len(R),len(D)
    CB = [[gal(z,-1) for z in row] for row in B]
    M = [row[1:] for row in CB]
    rhs = [[-add(C[k][j]*gal(P[i][j],-1) for j in range(s)) for k in range(n)]+[-CB[i][0]] for i in range(r)]
    sol = solve(M, rhs, K)
    b = [one] + [row[-1] for row in sol]
    base = [[zero]+[row[k] for row in sol] for k in range(n)]
    CC = [[gal(C[k][j],-1) for k in range(n)] for j in range(s)]
    rhs = [[-add(gal(P[i][j],-1)*B[i][0] for i in range(r)),-CC[j][0]] for j in range(s)]
    sol2 = solve([row[1:] for row in CC], rhs, K)
    f = [zero]+[row[0] for row in sol2]
    g = [one]+[row[1] for row in sol2]
    U = [[base[k][j]+f[k]*b[j] for j in range(l)] for k in range(n)]
    V = [[g[k]*b[j] for j in range(l)] for k in range(n)]
    eqs = []
    for i in range(n):
        for j in range(l):
            u,v = U[i][j],V[i][j]
            if v != zero:
                eqs.append((u,gal(v,2)/v-v,(gal(u,2)-u*u)/v))
    if all(d==zero for _,d,_ in eqs):
        t = next(z/u for u,d,z in eqs if u!=zero)/K.convert(2)
    else:
        a,b,c = next(q for q in eqs if q[1]!=zero)
        d,e,f = next(q for q in eqs if a*q[1]!=b*q[0])
        t = (c*e-b*f)/(a*e-b*d)/K.convert(2)
    assert all(U[i][j]+t*V[i][j]==H(row,col) for i,row in enumerate(R) for j,col in enumerate(D))
    print(f'p={p}, {r}x{s} Fourier overlap, recovered all {n*l} entries exactly ({monotonic()-start:.1f}s)',flush=True)

run(5, [0,2], [0,3])
run(7, [0,2,4], [0,2,4])
run(7, [1,5], [0,2,3,6])
run(11,[0,2,4,6,8],[0,1,3,4,8])
