"""Exact first classification of profile feasibility with an 11-cycle quotient.

Small-root cycle values are matched in mu_6, preserving multiplicities.
The 11-component of the phase at a cycle is a vertex label in Z/11.
An edge from a left cycle to a right cycle has label right-left. These 11
labels must be distinct. Parallel edges are therefore impossible.
"""
from itertools import combinations_with_replacement, permutations, product
from collections import Counter, defaultdict
from pathlib import Path
import json, time

ROOT=Path(__file__).resolve().parent
P={'E':(3,3,3,2),'F':(3,2,2,2,2)}

def configs(shape):
    groups=[]
    for p in sorted(set(shape),reverse=True):
        opts=[tuple((a+6*k//p)%6 for k in range(p)) for a in range(6//p)]
        groups.append(list(combinations_with_replacement(opts,shape.count(p))))
    return [tuple(c for g in parts for c in g) for parts in product(*groups)]

def multiplicities(config):
    return tuple(sum(r in c for c in config) for r in range(6))

def shifted(config,t):
    groups=defaultdict(list)
    for c in config: groups[len(c)].append(tuple(sorted((r+t)%6 for r in c)))
    return tuple(c for p in sorted(groups,reverse=True) for c in sorted(groups[p]))

def normalized_pair(a,b):
    return min((shifted(a,t),shifted(b,t)) for t in range(6))

def column_perms(shape):
    groups=[[i for i,s in enumerate(shape) if s==p] for p in sorted(set(shape),reverse=True)]
    return [tuple(i for group in choice for i in group)
            for choice in product(*(list(permutations(g)) for g in groups))]

def canonical_graph(edges, left_shape, right_shape):
    best=None; best_edges=None
    for columns in column_perms(right_shape):
        ci={old:new for new,old in enumerate(columns)}
        masks=[sum(1<<ci[k] for i,k,r in edges if i==j) for j in range(len(left_shape))]
        order=[]
        for p in sorted(set(left_shape),reverse=True):
            order += sorted((i for i,s in enumerate(left_shape) if s==p),key=lambda i:masks[i])
        key=tuple(masks[i] for i in order)
        if best is None or key<best:
            ri={old:new for new,old in enumerate(order)}
            best=key
            best_edges=tuple(sorted((ri[i],ci[k],r) for i,k,r in edges))
    return best,best_edges

def graphs(a,b):
    normalized=set()
    for ca in configs(P[a]):
        for cb in configs(P[b]):
            if multiplicities(ca)==multiplicities(cb):
                normalized.add(normalized_pair(ca,cb))
    out={}; matches=0; distinct_colored=set()
    for ca,cb in sorted(normalized):
        # Empty values have a unique empty matching.
        by_root=[]
        for r in range(6):
            L=[i for i,c in enumerate(ca) if r in c]
            R=[i for i,c in enumerate(cb) if r in c]
            by_root.append([tuple((i,k,r) for i,k in zip(L,rp)) for rp in permutations(R)])
        def rec(r,edges,occupied):
            nonlocal matches
            if r==6:
                matches+=1
                key,ce=canonical_graph(edges,P[a],P[b])
                distinct_colored.add(ce)
                record=out.setdefault(key,{'edges':ce,'left_shape':P[a],'right_shape':P[b],'colorings':set()})
                record['colorings'].add(ce)
                return
            for block in by_root[r]:
                pairs={(i,k) for i,k,v in block}
                if pairs.isdisjoint(occupied): rec(r+1,edges+block,occupied|pairs)
        rec(0,(),set())
    for record in out.values(): record['colorings']=sorted(record['colorings'])
    return out,{'phase_configurations_up_to_global_mu6_rotation':len(normalized),
                'simple_colored_matchings':matches,'canonical_uncolored_graphs':len(out),
                'canonical_colored_representations':len(distinct_colored)}

def find_labels(rec):
    ls,rs=rec['left_shape'],rec['right_shape']; L=len(ls); R=len(rs); V=L+R
    edges=[(i,L+k) for i,k,r in rec['edges']]
    adj=[[] for _ in range(V)]
    for i,k in edges: adj[i].append(k);adj[k].append(i)
    # Independent translations of the two sides and common multiplication
    # by a nonzero field element allow b_i=0,b_j=1,c_k=0 at any right vertex.
    k=L
    i,j=adj[k][:2]
    values=[None]*V; values[i]=0;values[j]=1;values[k]=0
    nodes=0
    def mask_for_current():
        mask=0
        for x,y in edges:
            if values[x] is not None and values[y] is not None:
                d=(values[y]-values[x])%11
                if mask>>d&1:return None
                mask|=1<<d
        return mask
    def search(mask):
        nonlocal nodes
        nodes+=1
        unset=[x for x in range(V) if values[x] is None]
        if not unset: return values[:]
        x=max(unset,key=lambda x:(sum(values[y] is not None for y in adj[x]),len(adj[x])))
        domain=range(11)
        if len(unset)==1:
            total=sum((len(adj[y]) if y>=L else -len(adj[y]))*values[y]
                      for y in range(V) if y!=x)
            coeff=len(adj[x]) if x>=L else -len(adj[x])
            domain=[(-total*pow(coeff%11,-1,11))%11]
        for value in domain:
            m=mask
            for y in adj[x]:
                if values[y] is None: continue
                d=(value-values[y])%11 if x>=L else (values[y]-value)%11
                if m>>d&1:break
                m|=1<<d
            else:
                values[x]=value
                result=search(m)
                if result is not None:return result
                values[x]=None
        return None
    mask=mask_for_current()
    labels=search(mask) if mask is not None else None
    return labels,nodes

def witness(rec,labels):
    L=len(rec['left_shape'])
    columns=[]
    for i,k,r in rec['edges']:
        b,c=labels[i],labels[L+k]
        j=(c-b)%11
        # Values are mu_11^b * mu_6^r and mu_11^c * mu_6^r.
        # H=(1, u, w), with u_j=mu_11^j and w_j=mu_11^c*mu_6^r;
        # then w/u is the left cycle vector.
        columns.append((j,(6*c+11*r)%66))
    columns.sort()
    w=[x[1] for x in columns]; w0=w[0]
    return {'modulus':66,'row_exponents':[[0]*11,[6*j for j in range(11)],[(x-w0)%66 for x in w]],
            'cycle_edges':rec['edges'],'vertex_labels':labels}

if __name__=='__main__':
    report={}
    for a,b in [('E','E'),('E','F'),('F','F')]:
        t=time.time(); gs,stats=graphs(a,b)
        print(a+b,stats,flush=True)
        records=[]
        for n,(key,rec) in enumerate(sorted(gs.items())):
            labels,nodes=find_labels(rec)
            record={'graph':key,**rec,'label_search_nodes':nodes,'feasible':labels is not None}
            if labels is not None:record['witness']=witness(rec,labels)
            records.append(record)
            print(' ',n+1,key,'feasible',labels is not None,'nodes',nodes,flush=True)
        stats['feasible_graphs']=sum(x['feasible'] for x in records)
        stats['elapsed_seconds']=time.time()-t
        report[a+b]={'summary':stats,'graphs':records}
        (ROOT/'results.json').write_text(json.dumps(report,indent=2)+'\n')
        print('FINISHED',a+b,stats,flush=True)
