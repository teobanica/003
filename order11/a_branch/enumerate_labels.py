"""Exhaust the finite modulo-11 labels on each of the six support graphs.

Counts are normalized label assignments, not Hadamard equivalence classes.
Fix b_i=0,b_j=1,c_k=0 at the first right vertex k and its first two neighbors.
For every other graph component, fix its first vertex label to zero.
"""
from pathlib import Path
from collections import Counter
import json,time
ROOT=Path(__file__).resolve().parent

def enumerate_labels(rec):
    L=len(rec['left_shape']); V=L+len(rec['right_shape'])
    edges=[(i,L+k) for i,k,r in rec['edges']]
    adj=[[] for _ in range(V)]
    for i,k in edges:adj[i].append(k);adj[k].append(i)
    remaining=set(range(V)); comps=[]
    while remaining:
        seen=set(); stack=[min(remaining)]
        while stack:
            x=stack.pop()
            if x in seen:continue
            seen.add(x); stack.extend(adj[x])
        comps.append(sorted(seen));remaining-=seen
    values=[None]*V
    k=L;i,j=adj[k][:2]
    values[i]=0;values[j]=1;values[k]=0
    for comp in comps:
        if not any(values[x] is not None for x in comp):values[comp[0]]=0
    fixed={x:values[x] for x in range(V) if values[x] is not None}
    mask=0
    for x,y in edges:
        if values[x] is not None and values[y] is not None:
            d=(values[y]-values[x])%11
            if mask>>d&1:raise ValueError('inconsistent normalization')
            mask|=1<<d
    solutions=[];nodes=0
    def search(mask):
        nonlocal nodes
        nodes+=1
        unset=[x for x in range(V) if values[x] is None]
        if not unset:
            solutions.append(values[:]);return
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
                if values[y] is None:continue
                d=(value-values[y])%11 if x>=L else (values[y]-value)%11
                if m>>d&1:break
                m|=1<<d
            else:
                values[x]=value;search(m);values[x]=None
    search(mask)
    return {'components':comps,'fixed_vertex_labels':fixed,'normalised_label_count':len(solutions),
            'search_nodes':nodes,'label_solutions':solutions}

if __name__=='__main__':
    source=json.loads((ROOT/'results.json').read_text());out={}
    for kind,data in source.items():
        out[kind]=[]
        for rec in data['graphs']:
            start=time.time();r=enumerate_labels(rec);r['graph']=rec['graph']
            out[kind].append(r)
            print(kind,rec['graph'],'components',len(r['components']),'labels',r['normalised_label_count'],
                  'nodes',r['search_nodes'],'seconds',round(time.time()-start,3),flush=True)
            (ROOT/'label_enumeration.json').write_text(json.dumps(out,indent=2)+'\n')
