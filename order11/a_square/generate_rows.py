"""All dephased mu_11 rows orthogonal to 1 and (eta^j)."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent
rows=[]
def visit(j,values,differences,row):
    if j==11:rows.append(tuple(row));return
    for y in range(1,11):
        if values>>y&1:continue
        d=(y-j)%11
        if differences>>d&1:continue
        visit(j+1,values|1<<y,differences|1<<d,row+[y])
visit(1,1,1,[0])
assert len(rows)==3441 and len(set(rows))==3441
for row in rows:
    assert sorted(row)==list(range(11))
    assert sorted((row[j]-j)%11 for j in range(11))==list(range(11))
text=str(len(rows))+'\n'+'\n'.join(' '.join(map(str,row)) for row in rows)+'\n'
(ROOT/'rows.txt').write_text(text)
report={'candidate_rows':len(rows),'sha256_rows_text':hashlib.sha256(text.encode()).hexdigest(),
        'fourier_candidate_indices_zero_based':[rows.index(tuple(k*j%11 for j in range(11))) for k in range(2,11)]}
(ROOT/'row_generation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
