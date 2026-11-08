from sympy import Poly, symbols, gcd, factorint
from math import comb
x=symbols('x')
def data(p,g):
 m=(p-1)//2
 q={pow(g,i,p) for i in range(m)}
 assert len(q)==m
 def chi(t):
  if t%p==0: return 0
  return 1 if t%p in q else -1
 c=[chi(pow(g,i,p)-1)-chi(pow(g,i,p)+1) for i in range(m)]
 c[0]-=1
 f=Poly(sum(b*x**j for j,b in enumerate(c)),x)
 assert f.eval(1)==0
 assert pow(g,int(f.diff().eval(1)),p)==4%p
 return m,f,c
for p,g in [(67,4),(71,49),(79,9),(103,25)]:
 m,poly,c=data(p,g)
 v=poly.exquo(Poly(x-1,x))
 for l in factorint(m):
  f=Poly(x**m-1,x,modulus=l)
  z=gcd(v.set_modulus(l),f)
  assert z.degree()==0
 print(p,'gcd one at',list(factorint(m)))
m,poly,c=data(43,9)
for l,d in [(3,2),(7,1)]:
 z=gcd(poly.set_modulus(l),Poly(x**21-1,x,modulus=l))
 assert z==Poly((x-1)**d,x,modulus=l)
q={pow(9,i,43) for i in range(21)}
assert [(x,(x+1)%43) for x in sorted(q) if (x+1)%43 in q and (x+1)*pow(x,-1,43)%43 in (1,6,36)]==[(16,17)]
assert (pow(9,12,43),pow(9,19,43))==(16,17)
print(43,'gcd multiplicities and the unique marked ordered pair verified')

# The printed short remainders use the square-free period in each characteristic.
def check_remainder(polynomial, period, ell, expected):
    observed = polynomial.set_modulus(ell).rem(Poly(x**period - 1, x, modulus=ell))
    assert observed == Poly(expected, x, modulus=ell), (period, ell, observed)

for p, g, cases in [
    (67, 4, [(11,3,-x**10-x**6+x**5+x**4-x**3-x), (3,11,-3*x**2+x+3)]),
    (71,49, [(7,5,-2*x**5-2*x**3+2*x**2+x+2), (5,7,-3*x**4-2*x**3-3)]),
]:
    m, polynomial, _ = data(p, g)
    quotient = polynomial.exquo(Poly(x-1,x))
    for period, ell, expected in cases:
        check_remainder(quotient, period, ell, expected)
_, polynomial, _ = data(31, 9)
assert polynomial.diff().eval(1) == -6
_, polynomial, _ = data(43, 9)
check_remainder(polynomial, 7, 3, -x**5-x**3-x**2+x-1)
check_remainder(polynomial, 3, 7, 3*x**2+2*x+2)
print("All displayed odd-order polynomial remainders agree.")
