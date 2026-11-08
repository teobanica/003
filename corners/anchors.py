from math import log,floor
def anchors(p):
 q={x*x%p for x in range(1,p)}; R=sorted(q-{1});N=sorted(set(range(1,p))-q)
 yy={b:sum(1<<i for i,y in enumerate(N) if (y-b)%p in q) for b in R}
 nmask=(1<<len(N))-1
 live={x:nmask for x in R}
 chosen=[];sizes=[]
 while any(v for v in live.values()):
  scores=[]
  for b in R:
   if b in chosen:continue
   pos=yy[b];neg=nmask^pos
   scores.append((sum((v & (0 if x==b else neg if (x-b)%p in q else pos)).bit_count() for x,v in live.items()),b))
  remaining,b=min(scores); chosen.append(b);sizes.append(remaining)
  for x in R:
   live[x] &= 0 if x==b else nmask^yy[b] if (x-b)%p in q else yy[b]
 for x in R:
  for y in N:
   assert any(b==x or ((y-b)*(x-b))%p in q for b in chosen)
 m=len(q);rho=0.5+(p**.5)/(2*(m-1))
 return chosen,sizes, min(m-1,1+floor(log(m*(m-1))/-log(rho)))
if __name__ == '__main__':
    chosen, remaining, bound = anchors(503)
    assert chosen == [2,16,231,104,286,392,144,66,197,134,376,4]
    assert remaining[-1] == 0
    print('p=503, twelve selected rows:', chosen)
