"""Independent symbolic/Fraction replay of the parent's fourteen direction signs."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import sympy as s

base=Path(__file__).resolve().parent
cert=json.loads((base.parent/'postclip_knot_direction.json').read_text())
r,rho,u,v=s.symbols('r rho u v');a,b=s.symbols('a b');sigma=1-rho
chi_1=1-r-r/rho;chi_2=1-r-r/sigma;c=1-2*r
v0=r*(1+sigma)/(sigma*(1-r))
h0=1-r*(1+rho)/sigma;h1=rho*(1-r)/sigma
T0=h0-r+v0*(chi_1-r);T1=h1+v0*chi_2;ut=s.cancel(-T0/T1)
H=s.cancel(T1*sigma*sigma*(r-1))
N_R=sigma*sigma*chi_2+rho*(sigma-2*r)*u+sigma*(sigma*c+(rho*c-r)*u)*v
N_C=rho*rho*chi_1+rho*rho*c*u+sigma*(rho-2*r+rho*(chi_2-r)*u)*v
L=(chi_1+c*u)*N_R+r*u*N_C
B_mix=(chi_2+c*v)*N_C+r*v*N_R
kv=s.diff(L,v)
ex={'negative_Kv_at_ut':-kv.subs(u,ut),
    'negative_Kvu_at_ut':-s.diff(kv,u).subs(u,ut),
    'negative_Kvuu':-s.diff(kv,u,2)}
den={'negative_Kv_at_ut':rho*H**2,
     'negative_Kvu_at_ut':-rho*H,
     'negative_Kvuu':s.Integer(1)}
for label,uv in [('ut',ut),('one',s.Integer(1))]:
    ex['negative_KBv_at_v0_'+label]=-s.diff(B_mix,v).subs({v:v0,u:uv})
    ex['negative_KBvv_'+label]=-s.diff(B_mix,v,2).subs(u,uv)
    den['negative_KBv_at_v0_'+label]=rho*(r-1)*H if label=='ut' else (rho-1)*(r-1)
    den['negative_KBvv_'+label]=-rho*H if label=='ut' else s.Integer(1)
assert s.factor(T1-(H/(sigma*sigma*(r-1))))==0
assert s.diff(s.diff(B_mix,v),u,2)==0 and s.diff(s.diff(B_mix,v,2),u,2)==0
# T1>0 follows from chi_2>=0, or T1>=h1+chi_2/k=(1-r)(1-2r)/r>0.
# Consequently H<0; every denominator above is positive on the whole domain.
assert s.factor(h1+chi_2/(r/(1-r))-(1-r)*(1-2*r)/r)==0

def convert(power,da,db):
    intermediate={}
    for i in range(da+1):
        for j in range(db+1):
            intermediate[i,j]=sum(power.get((k,j),F(0))*F(comb(i,k),comb(da,k)) for k in range(i+1))
    return {(i,j):sum(intermediate[i,l]*F(comb(j,l),comb(db,l)) for l in range(j+1))
            for i in range(da+1) for j in range(db+1)}

r0=s.Rational(30953996,10**8);rr=r0+(s.Rational(1,3)-r0)*a
out=[]
for half,pp in [('small',rr+(s.Rational(1,2)-rr)*b),
                ('large',s.Rational(1,2)+(s.Rational(1,2)-rr)*b)]:
    for name,value in ex.items():
        ent=cert['results'][half+'_'+name]
        recorded=s.sympify(ent['positive_denominator_to_justify'],locals={'r':r,'rho':rho})
        assert s.factor(recorded-den[name])==0
        numerator=s.cancel(value*den[name])
        assert s.fraction(numerator)[1]==1
        po=s.Poly(s.expand(numerator.subs({r:rr,rho:pp})),a,b)
        degrees=list(po.degree_list());assert degrees==ent['degree']
        vals=convert({idx:F(cv) for idx,cv in po.terms()},*degrees)
        assert len(vals)==ent['count'] and min(vals.values())==F(ent['minimum'])
        assert all(cv>0 for cv in vals.values())
        out.append({'half':half,'name':name,'count':len(vals),'all_positive':True})
report={'success':True,'subdivisions':0,'independent_symbolic_derivation':True,
        'independent_fraction_basis_conversion':True,
        'positive_denominators_from_T1':True,'coefficient_count':sum(chi_1['count'] for chi_1 in out),
        'results':out}
(base/'postclip-knot-direction-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
