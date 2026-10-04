"""Independent full-mass cap-only both-high replay, allowing negative s_1/s_2."""
import json
from pathlib import Path
from fractions import Fraction
from itertools import product
from math import comb
import sympy as s

folder=Path(__file__).parent
saved=json.loads((folder.parent/'fullmass_both_high.json').read_text())
r,rho,u,v,a,b,w,z=s.symbols('r rho u v a b w z');sigma=1-rho;c=1-2*r
chi_1=1-r-r/rho;chi_2=1-r-r/sigma
h_1=1-r*(1+sigma)/rho+sigma*(1-r)*v/rho
h_2=1-r*(1+rho)/sigma+rho*(1-r)*u/sigma
S=h_1-r+u*(chi_2+chi_1*v-r);T=h_2-r+v*(chi_1+chi_2*u-r)
Sv=s.diff(S,v);Tv=s.diff(T,v)
assert s.factor(Sv.subs(u,0)-sigma*(1-r)/rho)==0
assert s.factor(Sv.subs(u,1)-c/rho)==0
assert s.factor(Tv.subs(u,0)-(chi_1-r))==0
assert s.factor(Tv.subs(u,1)-(chi_1+chi_2-r))==0
assert s.factor(chi_1-r-(1-2*r-r/rho))==0
assert s.factor(chi_1+chi_2-r-(2-3*r-r/(rho*sigma)))==0
r0=s.Rational('0.30953996')
# chi_1-r <= (1-4r+2r^2)/(1-r)<0 for rho<=1-r.
outer=1-4*r+2*r*r
assert outer.subs(r,r0)<0
assert s.diff(outer,r)==-4+4*r
assert 2-7*r0<0
# No lower bound on chi_1+yu or chi_2+xv entered these endpoint arguments.

vs=s.factor((r-S.subs(v,0))/Sv)
vs0=vs.subs(u,0);vs1=vs.subs(u,1)
N=s.factor(s.diff(T,u)*Sv-Tv*s.diff(S,u))
assert s.factor(N-((1-r)**2-(chi_2-r)*(chi_1-r)+u*(rho*(1-r)*chi_1/sigma-chi_2*(chi_2-r))+v*(sigma*(1-r)*chi_2/rho-chi_1*(chi_1-r))))==0
assert s.factor(vs1-(r*(sigma+1/sigma)/c-2*rho))==0
exprs={'jacobian':N,'endpoint_r_minus_T':r-T.subs({u:1,v:vs1}),
       'vs0':vs0,'vs1':vs1,'vs0_upper':s.Rational(7,3)-vs0,
       'vs1_upper':s.Rational(7,3)-vs1}
rr=r0+(s.Rational(1,3)-r0)*a
pp=rr+(1-2*rr)*b
tot=0;checked={}
for name,expr in exprs.items():
    cert=saved[name]
    denom=s.sympify(cert['positive_denominator'],locals={'r':r,'rho':rho})
    numerator=s.factor(expr*denom)
    assert s.denom(numerator)==1
    nn,dd=s.fraction(s.factor(numerator.subs({r:rr,rho:pp,u:w,v:s.Rational(7,3)*z})))
    assert dd==s.sympify(cert['substitution_denominator']) and dd>0
    variables=(a,b,w,z) if name=='jacobian' else (a,b)
    poly=s.Poly(s.expand(nn),*variables)
    assert poly.degree_list()==tuple(cert['original_degree'])
    powers={ks:Fraction(cv) for ks,cv in poly.terms()};degrees=tuple(cert['degree'])
    vals=[]
    for ids in product(*[range(n+1) for n in degrees]):
        val=Fraction(0)
        for ks,co in powers.items():
            if not all(k<=i for k,i in zip(ks,ids)):continue
            factor=Fraction(1)
            for i,k,n in zip(ids,ks,degrees):factor*=Fraction(comb(i,k),comb(n,k))
            val+=co*factor
        assert val>0
        vals.append(val)
    assert len(vals)==cert['count'] and min(vals)==Fraction(cert['minimum'])
    tot+=len(vals);checked[name]=len(vals)

# Record denominator signs by their exact factorizations.
assert s.factor(s.sympify(saved['jacobian']['positive_denominator'])-rho*rho*sigma*sigma)==0
assert s.factor(s.sympify(saved['endpoint_r_minus_T']['positive_denominator'])-rho*sigma*sigma*c)==0
for name,expected in [('vs0',sigma*(1-r)),('vs1',sigma*c),
                      ('vs0_upper',3*sigma*(1-r)),('vs1_upper',3*sigma*c)]:
    assert s.factor(s.sympify(saved[name]['positive_denominator'])-expected)==0
report={'raw_model_identities_verified':True,
 'S_v_strictly_positive_full_rho_u_0_to_1':True,
 'T_v_strictly_negative_full_rho_u_0_to_1':True,
 'cap_lower_bounds_used':False,'jacobian_identity_verified':True,
 'all_positive_coefficient_counts':checked,'total_positive_coefficients':tot,
 'all_stored_minima_matched':True,'all_response_denominators_positive':True,'no_subdivision':True,
 'application':'For full rho in [r,1-r], r>=.30953996, 0<=u<=1, response-capped S>=r forces response-capped T<r with linear s_1/s_2. This excludes T=r in either mass orientation even after s_1/s_2 becomes negative.',
 'result':'passed'}
(folder/'fullmass-both-high-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
