"""Independently derive the T=0 bound and replay sequential Bernstein transforms."""
import sympy as s
from fractions import Fraction as F
from math import comb
from itertools import product
from pathlib import Path
import json

base=Path(__file__).resolve().parent
cert=json.loads((base/'m-chart-certificate.json').read_text())
r,rho,v,U=s.symbols('r rho v U');sigma=1-rho
R,W,Z=s.symbols('R W Z')
c=s.Rational(cert['r_lower'])
chi_1=1-r-r/rho;chi_2=1-r-r/sigma
rawh=1-r*(1+rho)/sigma+rho*(1-r)/sigma*U
s_1=chi_1+chi_2*U;s_2=chi_2+chi_1*v
rawT=rawh-r+v*(s_1-r)
ut=s.solve(rawT,U)[0]
k=r/(1-r);vk=r*(1+sigma)/(sigma*(1-r))


def sequential_bernstein(power,degrees):
    vals=dict(power)
    for axis,n in enumerate(degrees):
        other=[i for i in range(3) if i!=axis]
        nxt={}
        for fixed in product(*(range(degrees[i]+1) for i in other)):
            for j in range(n+1):
                mon=[0,0,0]
                for ax,value in zip(other,fixed):mon[ax]=value
                total=F(0)
                for i in range(j+1):
                    mon[axis]=i
                    total+=vals.get(tuple(mon),F(0))*F(comb(j,i),comb(n,i))
                mon[axis]=j;nxt[tuple(mon)]=total
        vals=nxt
    return vals


summary=[]
for entry in cert['cubes']:
    zero=entry['chart']=='BR_zero'
    h_1=1 if zero else 1-(r*(1+sigma)-sigma*(1-r)*v)/rho
    rawS=h_1-r+U*(s_2-r)
    rawM=(1-r)*rawS-(1-2*r)*rawT-r*(1+v)*rawh
    # Multiply M(T=0) by positive dT/dU before clearing positive denominators.
    val=s.factor(rawM.subs(U,ut)*s.diff(rawT,U))
    vs=vk+(1/k-vk)*Z if zero else k+(vk-k)*Z
    num,den=s.fraction(s.factor(val.subs(v,vs)))
    expected=rho*rho*(rho-1)**4*(r-1)**2*(r*r if zero else 1)
    assert s.factor(den-expected)==0
    poly=s.Poly(s.expand(num.subs({r:c+(s.Rational(1,3)-c)*R,rho:c+(s.Rational(1,2)-c)*W})),R,W,Z)
    actual={mon:F(cv) for mon,cv in poly.terms()}
    recorded={tuple(mon):F(cv) for mon,cv in entry['power_coefficients']}
    assert actual==recorded
    assert list(poly.degree_list())==entry['degrees']
    bern=sequential_bernstein(recorded,entry['degrees'])
    assert len(bern)==entry['bernstein_count']
    assert max(bern.values())==F(entry['maximum_bernstein_coefficient'])
    assert all(co<0 for co in bern.values())
    summary.append({'chart':entry['chart'],'count':len(bern),'strictly_negative':True,'max':str(max(bern.values()))})

# The elementary kappa_2=0 bound uses only r>=309/1000 and s_2<=6/25.
rr=F(309,1000)
bc0=1-F(32,25)*rr-2*rr*rr-rr/(1-rr)
assert bc0==F(-11629531,345500000)<0
report={'success':True,'independent_polynomial_derivation':True,
        'independent_sequential_bernstein_replay':True,'subdivisions':0,
        'cubes':summary,'BC_zero_negative_upper':str(bc0)}
(base/'m-chart-certificate-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
