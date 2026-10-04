"""Exact Fraction Krawczyk check of ONE tiny local stationary-root box.

This certifies a locally unique algebraic root, not a global scalar optimum.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
import sympy as sp
OUT=Path(__file__).resolve().parent
BASE=OUT.parent.parent
sys.path.insert(0,str(BASE/'baseline_audit'))
from six_lp_shared_certificate import I

r,rho,eta=sp.symbols('r rho eta');syms=[r,rho,eta]
meta=json.loads((OUT/'kink-polynomials.json').read_text())
eq=[sp.sympify(meta[k]) for k in ['N_1','N_2','H_root']]
root=json.loads((OUT/'kink-kkt.json').read_text())['variables']
center=[F(root[k][:52]) for k in ['r','rho','eta']]
radius=F(1,10**30)
box=[I(c-radius,c+radius) for c in center]

def peval(poly,vals):
    terms=sp.Poly(poly,*syms).terms()
    degrees=[max(rho[k] for rho,c in terms) for k in range(3)]
    powers=[]
    for val,d in zip(vals,degrees):
        ps=[I.co(1) if isinstance(val,I) else F(1)]
        for _ in range(d):ps.append(ps[-1]*val)
        powers.append(ps)
    out=I.co(0) if isinstance(vals[0],I) else F(0)
    for mon,coef in terms:
        term=F(coef)
        for k,rho in enumerate(mon):term=term*powers[k][rho]
        out=out+term
    return out

jpoly=[[sp.diff(f,x) for x in syms] for f in eq]
jcenter=sp.Matrix([[sp.Rational(peval(f,center)) for f in row] for row in jpoly])
cinv=jcenter.evalf(85).inv()
C=[[F(str(sp.N(cinv[i,j],55))) for j in range(3)] for i in range(3)]
J=[[peval(f,box) for f in row] for row in jpoly]
fc=[peval(f,center) for f in eq]
Y=[[I.co(1 if i==j else 0)-sum((C[i][k]*J[k][j] for k in range(3)),I.co(0))
    for j in range(3)] for i in range(3)]
contraction_norm=max(sum(max(abs(y.lo),abs(y.hi)) for y in row) for row in Y)
K=[]
for i in range(3):
    shift=-sum(C[i][j]*fc[j] for j in range(3))
    K.append(I.co(shift)+sum((Y[i][j]*I(-radius,radius) for j in range(3)),I.co(0)))
included=all(-radius<k.lo<=k.hi<radius for k in K)
contraction=contraction_norm<1
assert included and contraction
# Coarse rational bounds keep the saved certificate small and replayable.
qbound=F(1,10**20)
krbound=F(1,10**19)
assert contraction_norm<qbound
assert all(max(abs(k.lo),abs(k.hi))/radius<krbound for k in K)

report={'scope':'Exact local unique root of N_1=N_2=H_root=0 inside one tiny box; no global optimality claim',
        'center':[str(c) for c in center],'radius':str(radius),
        'preconditioner':[[str(c) for c in row] for row in C],
        'krawczyk_strict_inclusion':included,'infinity_norm_contraction':contraction,
        'contraction_rational_upper':str(qbound),'relative_K_offset_rational_upper':str(krbound),
        'floating_contraction':float(contraction_norm),'floating_max_relative_K_offset':float(max(max(abs(k.lo),abs(k.hi)) for k in K)/radius)}
(OUT/'local-root-certificate.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='preconditioner'},indent=2),flush=True)
