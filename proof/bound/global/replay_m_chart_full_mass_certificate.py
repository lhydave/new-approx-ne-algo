"""Independent polynomial derivation plus exact sequential coefficient replay."""
import sympy as s
from pathlib import Path
from fractions import Fraction as F
from math import comb
from itertools import product
import json

base=Path(__file__).resolve().parent
cert=json.loads((base/'m-chart-full-mass-certificate.json').read_text())
r,rho,v,U=s.symbols('r rho v U');sigma=1-rho
R,W,Z=s.symbols('R W Z');c=s.Rational(cert['r_lower'])
s_1=1-r-r/rho+(1-r-r/sigma)*U
s_2=1-r-r/sigma+(1-r-r/rho)*v
h_2=1-r*(1+rho)/sigma+rho*(1-r)/sigma*U
T=h_2-r+v*(s_1-r)
ut=s.solve(T,U)[0]
k=r/(1-r);vk=r*(1+sigma)/(sigma*(1-r));uC=r*(1+rho)/(rho*(1-r))


def sequential(power,degrees):
    out=dict(power)
    for axis,n in enumerate(degrees):
        other=[i for i in range(3) if i!=axis];nxt={}
        for fixed in product(*(range(degrees[i]+1) for i in other)):
            for j in range(n+1):
                mon=[0,0,0]
                for ax,val in zip(other,fixed):mon[ax]=val
                total=F(0)
                for i in range(j+1):
                    mon[axis]=i
                    total+=out.get(tuple(mon),F(0))*F(comb(j,i),comb(n,i))
                mon[axis]=j;nxt[tuple(mon)]=total
        out=nxt
    return out


summary=[]
for entry in cert['polynomials']:
    name=entry['name'];big=name.startswith('big_')
    if big:
        zero=name.endswith('zero')
        h_1=1 if zero else 1-(r*(1+sigma)-sigma*(1-r)*v)/rho
        S=h_1-r+U*(s_2-r)
        M=(1-r)*S-(1-2*r)*T-r*(1+v)*h_2
        expr=s.factor(M.subs(U,ut)*s.diff(T,U))
        vs=vk+(1/k-vk)*Z if zero else k+(vk-k)*Z
        expr=expr.subs(v,vs)
        mapping=s.Rational(1,2)+(s.Rational(1,2)-r)*W
        expected=rho*rho*(rho-1)**4*(r-1)**2*(r*r if zero else 1)
    else:
        if name.endswith('low'):vv,h_1=k,1-r/rho
        elif name.endswith('kink'):vv,h_1=vk,s.Integer(1)
        else:vv,h_1=1/k,s.Integer(1)
        expr=(1-r)*(h_1-r+uC*(s_2.subs(v,vv)-r))-r*(1+vv)
        mapping=r+(1-2*r)*W
        expected=-rho*rho*(rho-1) if name.endswith('high') else rho*rho*(rho-1)*(r-1)
    num,den=s.fraction(s.factor(expr))
    if name.endswith('high'):num=-num;den=-den
    assert s.factor(den-expected)==0
    poly=s.Poly(s.expand(num.subs(rho,mapping).subs(r,c+(s.Rational(1,3)-c)*R)),R,W,Z)
    power={idx:F(val) for idx,val in poly.terms()}
    saved={tuple(idx):F(val) for idx,val in entry['power_coefficients']}
    assert power==saved and list(poly.degree_list())==entry['degrees']
    bern=sequential(power,entry['degrees']);assert len(bern)==entry['bernstein_count']
    special=entry['special_negative_quadratic']
    positive=[idx for idx,val in bern.items() if val>=0]
    if special is None:assert not positive
    else:
        i=special['R_index'];j=special['Z_index']
        assert positive==[(i,1,j)]
        b0,b1,b2=[bern[(i,z,j)] for z in range(3)]
        assert [str(b0),str(6*b1),str(15*b2)]==special['quadratic_coefficients']
        margin=60*b0*b2-36*b1*b1
        assert margin==F(special['negative_discriminant_margin'])>0
        assert b0<0 and b2<0
        assert all(bern[(i,z,j)]<0 for z in range(3,7))
    summary.append({'name':name,'coefficients':len(bern),'negative_polynomial_certified':True,'special_quadratic':special is not None})

rr=F(cert['r_lower']);f=5*rr**4+2*rr**2-4*rr+1
assert f==F(cert['uniform_s_2_quartic_at_lower'])<0
report={'success':True,'independent_polynomial_derivation':True,'independent_sequential_replay':True,
        'subdivisions':0,'polynomials':summary,'uniform_s_2_quartic_negative':True}
(base/'m-chart-full-mass-certificate-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
