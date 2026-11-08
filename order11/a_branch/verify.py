"""Independent degree-matrix census and exact cyclotomic witness checks."""
from itertools import combinations,permutations,product
from pathlib import Path
from collections import defaultdict
import json
from functools import lru_cache
ROOT=Path(__file__).resolve().parent

def degree_graphs(ls,rs):
    R=len(rs);out=set();raw=0
    col_groups=[[j for j,d in enumerate(rs) if d==p] for p in sorted(set(rs),reverse=True)]
    cperms=[sum(x,()) for x in product(*(tuple(permutations(g)) for g in col_groups))]
    def canonical(rows):
        keys=[]
        for perm in cperms:
            transformed=[sum(((row>>old)&1)<<new for new,old in enumerate(perm)) for row in rows]
            key=tuple(v for d in sorted(set(ls),reverse=True)
                      for v in sorted(transformed[i] for i,p in enumerate(ls) if p==d))
            keys.append(key)
        return min(keys)
    def rec(i,capacities,rows):
        nonlocal raw
        if i==len(ls):
            if not any(capacities):raw+=1;out.add(canonical(rows))
            return
        for subset in combinations([j for j,c in enumerate(capacities) if c],ls[i]):
            mask=sum(1<<j for j in subset)
            if i and ls[i]==ls[i-1] and mask<rows[-1]:continue
            new=list(capacities)
            for j in subset:new[j]-=1
            if max(new,default=0)>len(ls)-i-1:continue
            rec(i+1,new,rows+[mask])
    rec(0,list(rs),[])
    return out,raw

def root_coloring(key,ls,rs):
    L=len(ls);shapes=ls+rs
    edges=[(i,L+k) for i,row in enumerate(key) for k in range(len(rs)) if row>>k&1]
    cycle_masks=[[sum(1<<((a+6*t//p)%6) for t in range(p)) for a in range(6//p)] for p in shapes]
    used=[0]*len(shapes);colors=[None]*len(edges)
    def search():
        candidate=None
        for e,(i,j) in enumerate(edges):
            if colors[e] is not None:continue
            opts=[]
            for r in ([0] if all(c is None for c in colors) else range(6)):
                bit=1<<r
                if used[i]&bit or used[j]&bit:continue
                if all(any(((used[x]|bit)&cm)==(used[x]|bit) for cm in cycle_masks[x]) for x in (i,j)):
                    opts.append(r)
            if not opts:return None
            if candidate is None or len(opts)<len(candidate[1]):candidate=(e,opts)
        if candidate is None:return colors[:]
        e,opts=candidate;i,j=edges[e]
        for r in opts:
            used[i]|=1<<r;used[j]|=1<<r;colors[e]=r
            ans=search()
            if ans is not None:return ans
            colors[e]=None;used[i]^=1<<r;used[j]^=1<<r
        return None
    return search()

def cycle_check(colors,shape,side):
    for j,p in enumerate(shape):
        vals=[r for i,k,r in colors if (i if side==0 else k)==j]
        possible=[{(a+6*t//p)%6 for t in range(p)} for a in range(6//p)]
        assert len(vals)==p and set(vals) in possible

def polynomial_division(a,b):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    quotient=[0]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and any(a):
        assert b[-1]==1
        k=len(a)-len(b);c=a[-1];quotient[k]=c
        for j,v in enumerate(b):a[k+j]-=c*v
        while len(a)>1 and a[-1]==0:a.pop()
    return quotient,a

@lru_cache(None)
def cyclotomic(n):
    poly=[-1]+[0]*(n-1)+[1]
    for d in range(1,n):
        if n%d==0:
            poly,rem=polynomial_division(poly,cyclotomic(d))
            assert not any(rem)
    return tuple(poly)

def check_family(rows,q_powers):
    phi=cyclotomic(66)
    for a,b in combinations(range(3),2):
        buckets=defaultdict(list)
        for j in range(11):buckets[q_powers[a][j]-q_powers[b][j]].append((rows[a][j]-rows[b][j])%66)
        for powers in buckets.values():
            poly=[0]*66
            for k in powers:poly[k]+=1
            assert not any(polynomial_division(poly,phi)[1])

def main():
    data=json.loads((ROOT/'results.json').read_text());labels=json.loads((ROOT/'label_enumeration.json').read_text())
    report={};brute_lines=[];all_labels=0;all_colors=0
    for kind,section in data.items():
        ls=section['graphs'][0]['left_shape'];rs=section['graphs'][0]['right_shape']
        graphs,raw=degree_graphs(ls,rs)
        passed={key for key in graphs if root_coloring(key,ls,rs) is not None}
        expected={tuple(rec['graph']) for rec in section['graphs']}
        assert passed==expected
        report[kind]={'degree_graph_isomorphism_classes':len(graphs),'phase_feasible_classes':len(passed),'agrees_with_cycle_matching_enumeration':True}
        print(kind,report[kind],flush=True)
        for rec,lab in zip(section['graphs'],labels[kind]):
            assert rec['graph']==lab['graph']
            for color in rec['colorings']:
                cycle_check(color,ls,0);cycle_check(color,rs,1);all_colors+=1
            L=len(ls);edges=[(i,L+k) for i,k,r in rec['edges']]
            for point in lab['label_solutions']:
                assert sorted((point[k]-point[i])%11 for i,k in edges)==list(range(11))
                all_labels+=1
            check_family(rec['witness']['row_exponents'],[[0]*11 for _ in range(3)])
            brute_lines.append(f"{kind} {len(ls)+len(rs)} {L} {len(lab['fixed_vertex_labels'])} {lab['normalised_label_count']}")
            brute_lines += [f"{i} {v}" for i,v in lab['fixed_vertex_labels'].items()]
            brute_lines += [f"{i} {k}" for i,k in edges]
            if len(lab['components'])==2:
                from compute import witness
                point=lab['label_solutions'][0];w=witness(rec,point)
                vertices=set(lab['components'][1])
                support=sorted((point[L+k]-point[i])%11 for i,k,r in rec['edges'] if i in vertices)
                qpow=[[0]*11,[0]*11,[int(j in support) for j in range(11)]]
                check_family(w['row_exponents'],qpow)
                cycles={}
                for side,shape in [(0,ls),(1,rs)]:
                    groups=[]
                    for z,p in enumerate(shape):
                        groups.append(sorted((point[L+k]-point[i])%11 for i,k,r in rec['edges'] if (i if side==0 else k)==z))
                    cycles['w_over_u' if side==0 else 'w']=groups
                family={'modulus':66,'row_exponents':w['row_exponents'],'q_exponents':qpow,'cycle_supports':cycles,'verified_by_exact_cyclotomic_remainders':True}
                (ROOT/'affine_family.json').write_text(json.dumps(family,indent=2)+'\n')
    report['verified_normalized_label_assignments']=all_labels
    report['verified_color_representations']=all_colors
    report['six_matrix_witnesses_and_symbolic_affine_family_exact']=True
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'brute_input.txt').write_text(str(sum(len(s['graphs']) for s in data.values()))+'\n'+'\n'.join(brute_lines)+'\n')
    print('Verified labels:',all_labels,'color representations:',all_colors,'and all cyclotomic identities.',flush=True)
if __name__=='__main__':main()
