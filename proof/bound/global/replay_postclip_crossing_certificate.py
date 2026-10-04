"""Independent C derivation and sequential exact Bernstein replay."""
import sympy as s
from fractions import Fraction as F
from math import comb
from itertools import product
from pathlib import Path
import json

base=Path(__file__).resolve().parent
cert=json.loads((base/'postclip-crossing-certificate.json').read_text())
r,rho,u=s.symbols('r rho u');R,W,Z=s.symbols('R W Z');sigma=1-rho
s_1=1-r-r/rho+(1-r-r/sigma)*u
h_2=1-r*(1+rho)/sigma+rho*(1-r)/sigma*u
expr=(r-s_1)*(2*r-s_1)-(h_2-r)*u*(1-r-r/rho)
k=r/(1-r);c=s.Rational(cert['r_lower'])


def sequential(vals,degrees):
    out=dict(vals)
    for axis,n in enumerate(degrees):
        other=[i for i in range(3) if i!=axis];nxt={}
        for fixed in product(*(range(degrees[i]+1) for i in other)):
            for j in range(n+1):
                mon=[0,0,0]
                for ax,val in zip(other,fixed):mon[ax]=val
                total=F(0)
                for i in range(j+1):
                    mon[axis]=i;total+=out.get(tuple(mon),F(0))*F(comb(j,i),comb(n,i))
                mon[axis]=j;nxt[tuple(mon)]=total
        out=nxt
    return out


summary=[]
for ent in cert['polynomials']:
    small=ent['half']=='small'
    pm=k+(s.Rational(1,2)-k)*W if small else s.Rational(1,2)+(s.Rational(1,2)-r)*W
    val=s.factor(expr.subs(u,k+(1-k)*Z).subs(rho,pm));num,den=s.fraction(val)
    expected=(3*W*r-W-2*r)**2*(3*W*r-W-4*r+2)**2 if small else (r-1)**2*(2*W*r-W-1)**2*(2*W*r-W+1)**2
    assert s.factor(den-expected)==0
    po=s.Poly(s.expand(num.subs(r,c+(s.Rational(1,3)-c)*R)),R,W,Z)
    power={idx:F(val) for idx,val in po.terms()}
    assert power=={tuple(idx):F(val) for idx,val in ent['power_coefficients']}
    assert list(po.degree_list())==ent['degree']
    bern=sequential(power,ent['degree'])
    assert len(bern)==ent['count'] and min(bern.values())==F(ent['minimum_bernstein'])
    assert all(co>0 for co in bern.values())
    summary.append({'half':ent['half'],'count':len(bern),'strictly_positive':True})
rr=F(cert['r_lower'])
assert 5*rr**4+2*rr**2-4*rr+1<0
assert F(319,6300)>0
# Independently clear the two rational derivative identities used by the barrier.
vv,ar,RR,XX,mm0,cc=s.symbols('v A RR XX m0 c')
nn0=r-ar;mm=mm0+XX*vv;nn=nn0+RR*vv;zz=nn/mm
ff=cc+vv*(-r-RR)-r*zz+vv*ar/zz
DD=RR*mm0-nn0*XX;CC=RR*(r+RR)-ar*XX
assert s.factor(s.diff(ff,vv)-(-CC/RR+DD*(ar*nn0/(RR*nn**2)-r/mm**2)))==0
assert s.factor(s.diff(ff,vv)-((r*zz-cc)/vv-(r+vv*ar/zz**2)*DD/mm**2)-ff/vv)==0
# The quartic is decreasing on the range: its derivative increases, but remains negative.
quartic=5*r**4+2*r**2-4*r+1
assert s.expand(s.diff(quartic,r,2)-(60*r**2+4))==0
assert s.diff(quartic,r).subs(r,s.Rational(1,3))==-s.Rational(52,27)
# The uniform zero-point budget has its minimum at the upper endpoint.
budget=(2*r-s.Rational(127,700))*(1-2*r)-r**2
assert budget.subs(r,s.Rational(1,3))==s.Rational(319,6300)
assert s.diff(budget,r).subs(r,s.Rational(3,10))<0
report={'success':True,'independent_polynomial_derivation':True,'independent_bernstein_replay':True,
        'subdivisions':0,'polynomials':summary,'quartic_strictly_negative':True,
        'two_derivative_identities':True,'uniform_zero_point_budget':True}
(base/'postclip-crossing-certificate-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
