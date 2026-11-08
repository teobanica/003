"""Integral dephased tangent matrix of a sixth-root Hadamard matrix."""
from pathlib import Path
from itertools import combinations
import numpy as np
home=Path(__file__).parent
E=np.loadtxt(home/'exponents.txt',dtype=int)
n=len(E)
pair=list(combinations(range(n),2))
u=np.array([[1,0],[0,1],[-1,1],[-1,0],[0,-1],[1,-1]],dtype=int)
def tangent(E):
    rows=[]
    for i,k in pair:
        weights=u[(E[i]-E[k])%6]
        assert np.all(weights.sum(axis=0)==0)
        for w in weights.T:
            R=np.zeros((n,n),int);R[i]=w;R[k]=-w
            rows.append(R[1:,1:].flatten())
    return np.array(rows)
