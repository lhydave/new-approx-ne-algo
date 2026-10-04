"""Independent short replay of root's T=0 and coarse-v stop lemmas."""
from pathlib import Path
from fractions import Fraction
from itertools import product
from math import comb
import sympy as s
import json

folder=Path(__file__).parent
saved=json.loads((folder.parent/'preclip_tzero_stop.json').read_text())
r,rho,u,v,a,b,w=s.symbols('r rho u v a b w');sigma=1-rho
chi_1=1-r-r/rho;chi_2=1-r-r/sigma
k=r/(1-r);v0=r*(1+sigma)/(sigma*(1-r))
h_2=1-r*(1+rho)/sigma+rho*(1-r)*u/sigma
h_1=1-r*(1+sigma)/rho+sigma*(1-r)*v/rho
s_1=chi_1+chi_2*u;s_2=chi_2+chi_1*v
S=h_1-r+u*(s_2-r);T=h_2-r+v*(s_1-r)
t1=s.diff(T,u);t0=T.subs(u,0)
ut=-t0/t1
Tzero=s.factor(2*r-S.subs(u,ut))
assert s.factor(T.subs(u,ut))==0
assert s.factor(Tzero-s.sympify(saved['expression'],locals={'r':r,'rho':rho,'v':v}))==0
den=s.sympify(saved['positive_denominator_to_justify'],locals={'r':r,'rho':rho,'v':v})
assert s.factor(den-rho*rho*sigma*sigma*t1)==0
# T1=rho(1-r)/sigma+v*chi_2>0 since chi_2>=1-3r>=0 canonically.
assert s.factor(t1-(rho*(1-r)/sigma+v*chi_2))==0
assert s.factor(chi_2-(1-3*r)-r*(2-1/sigma))==0

r0=s.Rational('0.30953996')
rr=r0+(s.Rational(1,3)-r0)*a
pp=rr+(s.Rational(1,2)-rr)*b
vv=(k+(v0-k)*w).subs({r:rr,rho:pp})
numer=s.factor(Tzero*den)
nn,dd=s.fraction(s.factor(numer.subs({r:rr,rho:pp,v:vv})))
assert s.factor(dd-s.sympify(saved['substitution_denominator'],locals={'a':a}))==0

def coefficients(poly,variables):
    degrees=poly.degree_list();powers={ks:Fraction(cv) for ks,cv in poly.terms()}
    out=[]
    for ids in product(*[range(n+1) for n in degrees]):
        val=Fraction(0)
        for ks,cv in powers.items():
            if not all(j<=i for j,i in zip(ks,ids)):continue
            f=Fraction(1)
            for i,j,n in zip(ids,ks,degrees):f*=Fraction(comb(i,j),comb(n,j))
            val+=cv*f
        assert val>0
        out.append(val)
    return degrees,out

poly=s.Poly(s.expand(nn),a,b,w)
deg,co=coefficients(poly,(a,b,w))
assert tuple(saved['degree'])==deg and saved['count']==len(co)==135
assert Fraction(saved['minimum'])==min(co)

vbound=s.Rational(4,5)
lower=s.factor(r-S.subs({u:k,v:vbound}))
ll=saved['canonical_lower_v']
assert s.factor(lower-s.sympify(ll['expression'],locals={'r':r,'rho':rho}))==0
ld=s.sympify(ll['positive_denominator'],locals={'r':r,'rho':rho})
assert s.factor(ld-5*rho*sigma*(1-r))==0
ln=s.factor(lower*ld)
pol=s.Poly(s.expand(ln.subs({r:rr,rho:pp})),a,b)
dg,lc=coefficients(pol,(a,b))
assert tuple(ll['degree'])==dg and ll['count']==len(lc)==15
assert Fraction(ll['minimum'])==min(lc)

# S_v>0 for u<=1: if chi_1<0, S_v>=sigma(1-r)/rho+chi_1=(1-2r)/rho.
assert s.factor(s.diff(S,v)-(sigma*(1-r)/rho+u*chi_1))==0
assert s.factor(sigma*(1-r)/rho+chi_1-(1-2*r)/rho)==0
margin=vbound+k-(1-2*r)/r
assert s.factor(s.diff(margin,r)-(1/(1-r)**2+1/r**2))==0
assert margin.subs(r,r0)>0
assert s.factor(margin.subs(r,r0)-s.Rational(ll['mass_margin_at_lower_r']))==0

# G=0 gives A=z*d, hence E=z[c-r(v+z)].
z,A,d,c=s.symbols('z A d c')
E=v*A+(c-v*(r+d))*z-r*z*z
assert s.expand(E.subs(A,z*d)-z*(c-r*(v+z)))==0
report={'T_zero_identity_verified':True,'T_zero_positive_denominator':'rho^2*sigma^2*T_u; T_u>0',
 'T_zero_all_135_coefficients_positive':True,
 'S_source_v_lower_bound':'v>4/5',
 'v_lower_bound_all_15_coefficients_positive':True,
 'S_monotonicity_identities_verified':True,
 'four_fifths_plus_k_exceeds_c_over_r':True,
 'source_budget_zero_envelope_identity_verified':True,
 'T_zero_stop':'T=0 forces S<2r; z<1 would force S>2r, contradiction.',
 'U_budget_zero_stop':'E=z[c-r(v+z)]<0 since v>4/5,z>=k and 4/5+k>c/r.',
 'result':'passed'}
(folder/'preclip-stop-events-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
