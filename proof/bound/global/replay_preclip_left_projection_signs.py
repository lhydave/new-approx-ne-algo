"""Independent Bernstein replay of the weak high-rho LEFT projection guards.

Reconstructs nine W/K signs, nine W/B_mix signs and five source-direction signs. No
production code is executed and no boxes are searched. N_C has no sign premise.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import sympy as s
import json,hashlib

folder=Path(__file__).resolve().parent.parent
r,rho,u,v,a,b=s.symbols('r rho u v a b');sigma=1-rho;c=1-2*r
k=r/(1-r);h=(1-r)/r;chi_1=1-r-r/rho;chi_2=1-r-r/sigma
h_1=1-r*(1+sigma)/rho+sigma*(1-r)*v/rho
h_2=1-r*(1+rho)/sigma+rho*(1-r)*u/sigma
s_1=chi_1+chi_2*u;s_2=chi_2+chi_1*v
S=h_1-r+(s_2-r)*u;T=h_2-r+(s_1-r)*v
N_R=sigma*sigma*chi_2+rho*(sigma-2*r)*u+sigma*sigma*c*v+rho*sigma*(chi_1-r)*u*v
N_C=rho*rho*chi_1+sigma*(rho-2*r)*v+rho*rho*c*u+rho*sigma*(chi_2-r)*u*v
K=s.expand((chi_1+c*u)*N_R+r*u*N_C)
W=s.factor((1-r)*S+r*T-r)
assert s.diff(K,v,2)==0 and s.diff(W,v,2)==0
K0=s.factor(K.subs(v,0));K1=s.factor(s.diff(K,v))
W0=s.factor(W.subs(v,0));W1=s.factor(s.diff(W,v))
NW=s.factor(W0*K1-W1*K0)
assert s.Poly(NW,u).degree()<=3
assert s.Poly(W1,u).degree()<=1 and s.Poly(K1,u).degree()<=2
local_coeff=[s.diff(NW,u,j).subs(u,1)/s.factorial(j)for j in range(4)]
left_signs={
    'negative_NWB'+str(j):-sum(local_coeff[t]*(h-1)**t*s.Rational(comb(j,t),comb(3,t))for t in range(j+1))
    for j in range(4)}
left_signs.update({'negative_Wv1':-W1.subs(u,1),'negative_Wvuh':-W1.subs(u,h),
    'K1_at_one':K1.subs(u,1),'K1_first':s.diff(K1,u).subs(u,1),
    'K1_second':s.diff(K1,u,2)/2})
left_denominators={
    'negative_NWB0':rho*rho*sigma,'negative_NWB1':3*rho*rho*r*sigma,
    'negative_NWB2':3*rho*rho*r*r*sigma,'negative_NWB3':rho*rho*r*r*sigma,
    'negative_Wv1':rho*sigma,'negative_Wvuh':rho*r*sigma,
    'K1_at_one':rho,'K1_first':rho,'K1_second':s.Integer(1)}

hc_u=rho*(1-r)/sigma;T_u=hc_u+chi_2*v
B=s.factor(T_u*(r-s_1)-(h_2-r)*(r-s_2))
D=s.factor(hc_u*(r-chi_1)+chi_2*(h_2.subs(u,0)-r))
assert s.Poly(B,u,v).degree(u)<=1 and s.Poly(B,u,v).degree(v)<=1
assert s.factor(s.diff((h_2-r)/(r-s_1),u)-D/(r-s_1)**2)==0
directions={}
denom_directions={}
for lu,uu in [('lo',k),('hi',1)]:
    for lv,vv in [('lo',k),('hi',h)]:
        directions['negative_B_'+lu+lv]=-B.subs({u:uu,v:vv})
denom_directions={'negative_B_lolo':sigma*sigma*(1-r)**2,
    'negative_B_lohi':sigma*sigma,'negative_B_hilo':sigma*sigma*(1-r),
    'negative_B_hihi':r*sigma*sigma,'negative_D':sigma*sigma}
directions['negative_D']=-D

# The opposite branch has no independent u=1 entry elsewhere: reconstruct
# the three shifted-B_mix coefficient signs in this same auxiliary verifier.
B_mix=s.expand((chi_2+c*v)*N_C+r*v*N_R)
assert s.Poly(B_mix,v).degree()<=2
below_coeff=[s.diff(B_mix,v,j).subs(v,0)/s.factorial(j)for j in range(3)]
below_numerators=[
    s.factor(below_coeff[0]*W1**2-below_coeff[1]*W0*W1+below_coeff[2]*W0**2),
    s.factor(below_coeff[1]*W1-2*below_coeff[2]*W0),s.factor(below_coeff[2])]
below_signs={}
for label,expr in zip(('N0','N1','N2'),below_numerators):
    degree=s.Poly(expr,u).degree()
    shifted=[s.diff(expr,u,j).subs(u,1)/s.factorial(j)for j in range(degree+1)]
    for j in range(degree+1):
        below_signs[label+'B'+str(j)]=sum(shifted[t]*(h-1)**t*s.Rational(comb(j,t),comb(degree,t))for t in range(j+1))
below_denominators={
    'N0B0':rho*rho*sigma**3,'N0B1':3*rho*rho*r*sigma**3,
    'N0B2':3*rho*rho*r*r*sigma**3,'N0B3':rho*rho*r*r*sigma**3,
    'N1B0':rho*sigma*sigma,'N1B1':2*rho*r*sigma*sigma,'N1B2':rho*r*r*sigma*sigma,
    'N2B0':s.Integer(1),'N2B1':r}

rr=s.Rational(619,2000)+(s.Rational(1,3)-s.Rational(619,2000))*a
checks=[]
def replay(relative,expressions,denominators,mapping):
    path=folder/relative;cert=json.loads(path.read_text())
    assert set(cert)==set(expressions)
    count=0
    for name,expr in expressions.items():
        saved=cert[name];den=denominators[name]
        # Each independent expected denominator is a product of positive
        # rho,sigma,r,1-r and positive integer constants on the stated domain.
        parsed=s.sympify(saved['positive_denominator'],locals={'rho':rho,'r':r})
        assert s.factor(parsed-den)==0
        numerator=s.factor(expr*den)
        assert s.denom(numerator)==1
        poly=s.Poly(s.expand(numerator.subs(mapping)),a,b)
        da,db=saved['degree'];assert da>=poly.degree(a)and db>=poly.degree(b)
        terms={ij:F(cv)for ij,cv in poly.terms()};values=[]
        for i in range(da+1):
            for j in range(db+1):
                z=sum(cv*F(comb(i,t),comb(da,t))*F(comb(j,w),comb(db,w))
                    for(t,w),cv in terms.items()if t<=i and w<=j)
                assert z==F(saved['coefficients'][i][j]) and z<0
                values.append(z);count+=1
        assert max(values)==F(saved['max'])
        checks.append({'name':name,'coefficients':len(values),'strictly_negative':True})
    return {'certificate':relative,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'coefficients_verified':count,'polynomials_verified':len(expressions)}

large_map={r:rr,rho:s.Rational(1,2)+(s.Rational(1,2)-rr)*b}
full_map={r:rr,rho:rr+(1-2*rr)*b}
results=[replay('kink/full_rho_w_signs.json',left_signs,left_denominators,large_map),
    replay('kink/full_rho_w_below_signs.json',below_signs,below_denominators,large_map),
    replay('kink/preclip_full_mass_direction_signs.json',directions,denom_directions,full_map)]
report={'success':True,'independent_symbolic_reconstruction':True,
    'independent_Fraction_Bernstein_replay':True,'positive_denominators_audited':True,
    'N_C_sign_assumed':False,'weak_LEFT_u1_barrier_verified':True,
    'weak_BELOW_u1_barrier_verified':True,
    'full_mass_fixed_v_budget_directions_verified':True,
    'scope':'r in [619/2000,1/3]; W/K barrier rho>=1/2; B,D full rho in [r,1-r]; u<=1 for directions',
    'certificates':results,'checks':checks}
out=Path(__file__).with_name('preclip-left-projection-signs-independent-replay.json')
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items()if k!='checks'},indent=2))
