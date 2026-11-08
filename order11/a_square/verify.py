"""Independent exhaustive certification of the all-11-cycle square completion.

Uses all 10! normalized permutations to verify the row list, then builds
compatibility by a direct set test. A 9-clique containing a degree-8 vertex
would require its entire neighborhood to be complete. Low degrees and one
missing neighbor edge per remaining non-Fourier vertex exclude every such
completion, without a general-purpose clique search.
"""
from pathlib import Path
from itertools import permutations,combinations
from collections import Counter
import json,time
ROOT=Path(__file__).resolve().parent

def main():
    t=time.time()
    rows=[tuple(map(int,s.split())) for s in (ROOT/'rows.txt').read_text().splitlines()[1:]]
    independent=[]
    for tail in permutations(range(1,11)):
        used={0}
        for j,y in enumerate(tail,1):
            d=(y-j)%11
            if d in used:break
            used.add(d)
        else:independent.append((0,)+tail)
    assert independent==rows and len(rows)==3441
    print('Independent 10! permutation enumeration agrees:',len(rows),flush=True)
    adj=[set() for _ in rows];edge_list=[]
    for i,x in enumerate(rows):
        for j in range(i):
            y=rows[j]
            if len({(a-b)%11 for a,b in zip(x,y)})==11:
                adj[i].add(j);adj[j].add(i);edge_list.append((j,i))
    fourier={rows.index(tuple(k*j%11 for j in range(11))) for k in range(2,11)}
    low=[];missing=[];partial_hist=Counter()
    for i in range(len(rows)):
        if i in fourier:continue
        neighbors=sorted(adj[i])
        if len(neighbors)<8:low.append(i)
        else:
            assert len(neighbors)==8
            absent=next(((a,b) for a,b in combinations(neighbors,2) if b not in adj[a]),None)
            assert absent is not None
            missing.append({'vertex':i,'nonadjacent_neighbors':absent})
        # Largest all-A partial extension containing this nonlinear row.
        size=0
        for r in range(len(neighbors),0,-1):
            if any(all(b in adj[a] for a,b in combinations(choice,2))
                   for choice in combinations(neighbors,r)):
                size=r;break
        partial_hist[3+size]+=1
    assert all(b in adj[a] for a,b in combinations(fourier,2))
    completed=json.loads((ROOT/'completions.json').read_text())
    assert len(completed['completion_cliques'])==1
    assert set(completed['completion_cliques'][0])==fourier
    assert completed['edges']==len(edge_list)
    report={'independently_generated_rows':len(rows),'compatibility_edges':len(edge_list),
            'degree_histogram':dict(sorted(Counter(map(len,adj)).items())),
            'nonfourier_vertices_excluded_by_degree':len(low),
            'nonfourier_degree8_vertices_with_missing_neighbor_edge':len(missing),
            'fourier_vertices_zero_based':sorted(fourier),
            'fourier_vertices_form_K9':True,'only_square_completion_is_fourier':True,
            'maximal_all_A_extension_size_by_nonfourier_row':dict(sorted(partial_hist.items())),
            'elapsed_seconds':round(time.time()-t,3)}
    certificate={'low_degree_vertices':low,'degree8_missing_edges':missing,'edges':edge_list}
    (ROOT/'independent_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'exclusion_certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
