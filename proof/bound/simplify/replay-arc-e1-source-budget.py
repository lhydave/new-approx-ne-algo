"""Independent replay of the one-cube E(1) source-budget certificate."""
import json
from pathlib import Path
from fractions import Fraction
from itertools import product
from math import comb
import sympy as s

folder=Path(__file__).parent
saved=json.loads((folder/'arc-e1-source-budget-certificate.json').read_text())
r,rho,u,v,a,b,t=s.symbols('r rho u v a b t');sigma=1-rho;c=1-2*r
chi_1=1-r-r/rho;chi_2=1-r-r/sigma
s_1=chi_1+chi_2*u;s_2=chi_2+chi_1*v
h_1=1-(r*(1+sigma)-sigma*(1-r)*v)/rho
h_2=1-(r*(1+rho)-rho*(1-r)*u)/sigma
S=h_1-r+u*(s_2-r);T=h_2-r+v*(s_1-r)
G=(S-r)*(h_2-r)-(r-T)*(r-s_1)
E1=1-3*r+v*(h_2+s_1-3*r)
N_R=sigma*sigma*chi_2+rho*(sigma-2*r)*u+sigma*sigma*c*v+rho*sigma*(chi_1-r)*u*v
N_C=rho*rho*chi_1+sigma*(rho-2*r)*v+rho*rho*c*u+rho*sigma*(chi_2-r)*u*v
L=(chi_1+c*u)*N_R+r*u*N_C
L_scaled_0=rho*L.subs(v,0);v_den=-rho*s.diff(L,v)
target=s.sympify(saved['target'],locals={'r':r,'rho':rho,'u':u})
assert s.factor(target-rho*rho*sigma*sigma*v_den*(-E1-2*G).subs(v,L_scaled_0/v_den))==0
assert s.Poly(target,r,rho,u).degree_list()==(4,7,4)
mapping={r:s.Rational(619,2000)+(s.Rational(1,3)-s.Rational(619,2000))*a,
         rho:s.Rational(3,10)+s.Rational(1,5)*b,
         u:s.Rational(2,5)+s.Rational(3,5)*t}
poly=s.Poly(s.expand(target.subs(mapping)),a,b,t)
ns=tuple(saved['bernstein_degree'])
power={ks:Fraction(cv) for ks,cv in poly.terms()}
co={}
for ids in product(*[range(n+1) for n in ns]):
    # Direct tensor formula, independent of generator's axis transformations.
    expected=Fraction(0)
    for ks,cv in power.items():
        if not all(k<=i for k,i in zip(ks,ids)):continue
        factor=Fraction(1)
        for i,k,n in zip(ids,ks,ns):factor*=Fraction(comb(i,k),comb(n,k))
        expected+=cv*factor
    key=','.join(map(str,ids))
    assert expected==Fraction(saved['coefficients'][key])
    co[ids]=expected
groups={tuple(g['indices']):g for g in saved['groups']}
positive_ungrouped=0
for i in range(ns[0]+1):
    for j in range(ns[1]+1):
        bb=[co[i,j,k] for k in range(5)]
        if(i,j) in groups:
            assert bb[0]>0 and bb[1]>0 and bb[2]>0 and bb[4]>0 and bb[3]<0
            assert 3*bb[2]*bb[4]-2*bb[3]**2>0
            assert str(3*bb[2]*bb[4]-2*bb[3]**2)==groups[i,j]['gap']
            positive_ungrouped+=2
        else:
            assert all(cv>0 for cv in bb)
            positive_ungrouped+=5
assert len(groups)==6 and positive_ungrouped==582
assert min(cv for (i,j,k),cv in co.items() if (i,j) not in groups or k<2)>Fraction(7,10**6)
assert all(Fraction(g['gap'])>Fraction(9,10**11) for g in groups.values())

# The final z<1 implication uses this elementary exact identity.
z,A,vv,rr,cc,dd=s.symbols('z A v r c d')
E=vv*A+(cc-vv*(rr+dd))*z-rr*z*z
Eone=E.subs(z,1)
assert s.expand(E-z*Eone-(1-z)*(vv*A+rr*z))==0
# Verify only the audited corrected derivative reduction, modulo E=0.
zz,mm,vv,rr,cc,ar,bb,yy,hh,kappa,GG,DD=s.symbols('z m v r c A b chi_2 h1 kappa G D')
LL=vv*bb-cc
Fprime=-bb-kappa*vv*yy-kappa*vv*hh/zz+ar/zz-(rr+vv*ar/zz**2)*(DD-zz*GG)/mm
C0=mm*vv*ar-mm*kappa*vv**2*hh+2*vv**2*ar*GG-vv*LL*DD
C1=-mm*vv*bb-mm*kappa*vv**2*yy-2*rr*vv*DD-vv*LL*GG
diff=s.factor((mm*vv*zz*Fprime-C0-C1*zz)*zz)
assert s.rem(s.Poly(diff,zz),s.Poly(rr*zz**2+LL*zz-vv*ar,zz)).as_expr()==0
report={'independent_model_identity_verified':True,
 'all_600_bernstein_coefficients_recomputed':True,
 '582_ungrouped_coefficients_strictly_positive':True,
 'six_positive_quadratic_discriminants':True,
 'one_fixed_cube_no_subdivision':True,
 'z_less_than_one_identity_verified':True,
 'corrected_zero_derivative_identity_verified':True,'result':'passed'}
(folder/'arc-e1-source-budget-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
