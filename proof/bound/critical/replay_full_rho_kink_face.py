"""Independent coefficient/identity replay for the full-rho oriented kink proof.

Imports no producer. Uses axis-by-axis Fraction Bernstein conversion and an
independent rational interval type for the algebraic-root quotient.
"""
from pathlib import Path
from fractions import Fraction as F
from dataclasses import dataclass
import itertools,json,math
import sympy as s
OUT=Path(__file__).resolve().parent
r,rho,eta=s.symbols('r rho eta');t_1,t_2,t_3=s.symbols('t_1 t_2 t_3')
m=json.loads((OUT/'kink-polynomials.json').read_text())
b=json.loads((OUT/'below-alignment-polynomials.json').read_text())
N_1,N_2,H_root=[s.sympify(m[k]) for k in ['N_1','N_2','H_root']]
S,T,s_1,s_2,A_0=[s.sympify(m['cap_expressions'][k]) for k in ['S','T','s_1','s_2','alpha']]
aa,bb,cc=[s.Poly(N_1,eta).nth(k) for k in [2,1,0]]
delta=s.cancel((bb*bb-4*aa*cc)/(rho-1)**2)
A_eta_num=s.cancel(s.fraction(s.factor(s.diff(A_0,eta)))[0]/(rho-1)**2)
N_B=s.sympify(b['N_B']);G_B=s.sympify(b['N_2_etaB_numerator_factors']['factors'][-1]['polynomial'])
R_resultant=s.sympify(json.loads((OUT/'kink-two-equation-resultant.json').read_text())['factors'][-1]['polynomial'])
lo=s.Rational(30953996,10**8)
checks=[]
def falling_ratio(k,j,d):
    z=F(1)
    for a in range(j):z*=F(k-a,d-a)
    return z
def check(label,expr,mapping,variables,degrees=None,weak=False):
    poly=s.Poly(s.expand(expr.subs(mapping,simultaneous=True)),*variables)
    deg=tuple(poly.degree_list()) if degrees is None else tuple(degrees)
    assert all(a<=d for a,d in zip(poly.degree_list(),deg))
    coeff={mon:F(a) for mon,a in poly.terms()}
    for axis,d in enumerate(deg):
        nxt={}
        for mon,a in coeff.items():
            j=mon[axis]
            for k in range(j,d+1):
                key=list(mon);key[axis]=k;key=tuple(key)
                nxt[key]=nxt.get(key,F(0))+a*falling_ratio(k,j,d)
        coeff=nxt
    vals=[coeff.get(idx,F(0)) for idx in itertools.product(*[range(d+1) for d in deg])]
    assert min(vals)>=0 if weak else min(vals)>0,(label,float(min(vals)))
    checks.append(dict(label=label,count=len(vals),weak=weak,numeric_min=float(min(vals))))
def num(expr):
    n,d=s.fraction(s.factor(expr));sg=s.sign(d.subs({r:s.Rational(31,100),rho:s.Rational(55,100),eta:s.Rational(6,10)}))
    assert sg in [-1,1]
    return s.expand(n*sg)

# Canonical branch, including the missing crossing semantics.
rx=lo+(s.Rational(1,3)-lo)*t_1;px=rx+(s.Rational(1,2)-rx)*t_2
T0=s.factor(s.limit(T,eta,s.oo));T1=s.factor((T-T0)*eta);etaT=T1/(r-T0)
sigma=1-rho;D_c=rho+(1-2*rho)*eta;p_dagger=rho*(1-eta)/D_c;q_dagger=r*(2-rho)/(rho+2*sigma*r)
J_C=(r*(rho+eta)-rho*eta)/D_c;J_R=(r*(3-2*rho)-sigma)/(rho+2*sigma*r)
base=(1-p_dagger)*(1-q_dagger)*r
Rk=base+q_dagger*(1-p_dagger*r/rho)-J_R;Ck=base+p_dagger*(1-q_dagger*r/sigma)-J_C
Rl=(1-p_dagger)*r;Cl=(1-p_dagger)*r+p_dagger-J_C
assert s.cancel((Cl-r)*(Rk-r)-(Rl-r)*(Ck-r)-N_1/s.sympify(m['D_1']))==0
xi=sigma*(1-r)/(sigma+r);v_0=(1-xi)/xi
h_2=1-r/sigma+(rho/sigma)*(1-eta-r)/eta;K_0=1-r/xi+v_0*s_1
Psi_A_0=A_0*(1-A_0)*K_0+(1-A_0)**2*v_0*h_2-r*A_0-r*(1-A_0)*v_0
# Independently reconstruct the defining cap and root-polynomial bridges.
# Stored polynomial coefficients are proof inputs, not assumed identities.
u=(1-eta)/eta
chi_1=1-r-r/rho;chi_2=1-r-r/sigma
sr_expected=chi_1+chi_2*u
sc_expected=chi_2+chi_1*v_0
S_expected=1-r+u*(sc_expected-r)
T_expected=h_2-r+v_0*(sr_expected-r)
A_expected=(r-T_expected)/(S_expected-T_expected)
for actual,expected in [(s_1,sr_expected),(s_2,sc_expected),(S,S_expected),(T,T_expected),(A_0,A_expected)]:
    assert s.cancel(actual-expected)==0
D_2=rho*eta*sigma*(1-r)*(rho*eta*sigma**2*(1-r)*(S-T))**2
assert s.cancel(Psi_A_0*D_2-N_2)==0
assert s.Poly(s.expand(H_root-(s.diff(N_1,rho)*s.diff(N_2,eta)-s.diff(N_1,eta)*s.diff(N_2,rho))),r,rho,eta).is_zero
res=s.resultant(N_1,N_2,eta)
assert s.Poly(s.expand(res-sigma**6*rho**2*r**2*R_resultant),r,rho).is_zero
for label,expr in [('T1>0',num(T1)),('r-T0>0',num(r-T0)),
 ('Ck numerator decreases eta',-s.diff(num(Ck-r),eta)),('Ck<r at T=r',num((r-Ck).subs(eta,etaT)))]:
    check(label,expr,{r:rx,rho:px},[t_1,t_2])
px=rx+(s.Rational(21,50)-rx)*t_2
check('rho<=.42 discriminant exclusion',-(bb*bb-4*aa*cc),{r:rx,rho:px},[t_1,t_2])
rx=s.Rational(8,25)+(s.Rational(1,3)-s.Rational(8,25))*t_1
px=rx+(s.Rational(1,2)-rx)*t_2;ex=rx+(1-2*rx)*t_3
check('canonical r>=.32 exclusion',-N_1,{r:rx,rho:px,eta:ex},[t_1,t_2,t_3])
rx=lo+(s.Rational(8,25)-lo)*t_1;px=s.Rational(21,50)+s.Rational(2,25)*t_2
for label,expr in [('a<0',-aa),('N_1(.5)<0',-N_1.subs(eta,s.Rational(1,2))),
 ('N_1_eta(.5)>0',s.diff(N_1,eta).subs(eta,s.Rational(1,2))),('delta_rho>0',s.diff(delta,rho)),
 ('delta_r<0',-s.diff(delta,r)),('canonical A_0_eta implication',A_eta_num-6*delta),
 ('S_eta>0',num(r-s_2)),('T_eta<0',num(-s.diff(T,eta)*eta**2)),
 ('S(.5)>r',num(S.subs(eta,s.Rational(1,2))-r)),('T(.5)<r',num(r-T.subs(eta,s.Rational(1,2)))),
 ('canonical R_resultant_r<0',-s.diff(R_resultant,r))]:check(label,expr,{r:rx,rho:px},[t_1,t_2])

# High-rho determinant face, with only previously checked broad domains.
px=s.Rational(1,2)+(1-lo-s.Rational(1,2))*t_2
for label,expr in [('high a<0',-aa),('high N_1(.5)<0',-N_1.subs(eta,s.Rational(1,2))),
 ('high N_1_eta(.5)>0',s.diff(N_1,eta).subs(eta,s.Rational(1,2))),('high delta_rho>0',s.diff(delta,rho)),
 ('high delta_r<0',-s.diff(delta,r)),('high A_0_eta implication',A_eta_num-2*N_1.subs(eta,1-r)),
 ('high N_1_eta(1-r)>0',s.diff(N_1,eta).subs(eta,1-r)),('high S_eta>0',num(r-s_2)),
 ('high T_eta<0',num(-s.diff(T,eta)*eta**2)),('high S(.5)>r',num(S.subs(eta,s.Rational(1,2))-r)),
 ('high T(.5)<r',num(r-T.subs(eta,s.Rational(1,2)))),('high R_resultant_r<0',-s.diff(R_resultant,r))]:
    check(label,expr,{r:rx,rho:px},[t_1,t_2])
rx=s.Rational(8,25)+(s.Rational(1,3)-s.Rational(8,25))*t_1
px=s.Rational(1,2)+(1-rx-s.Rational(1,2))*t_2
check('high r>=.32 eta monotonicity',s.diff(N_1,eta).subs(eta,1-r),{r:rx,rho:px},[t_1,t_2])
check('high r>=.32 endpoint exclusion',-N_1.subs(eta,1-r),{r:rx,rho:px},[t_1,t_2],[8,8])

# Opposite knot orientation: independently verify the rational substitution.
Fb=(q_dagger*(1-r)-J_R)*(Ck-r)+q_dagger*r*(Rk-r)
assert s.cancel(Fb-N_B/(sigma*D_c*(rho+2*sigma*r)**2))==0
np=s.Poly(N_B,eta);denB=np.nth(1);numB=-np.nth(0)
np2=s.Poly(N_2,eta)
NE=sum(np2.nth(j)*numB**j*denB**(3-j) for j in range(4))
assert s.Poly(s.expand(NE-sigma**3*r*G_B),r,rho).is_zero
rx=lo+(s.Rational(1,3)-lo)*t_1;px=s.Rational(1,2)+(1-rx-s.Rational(1,2))*t_2
for label,expr in [('N_B(.5)<0',-N_B.subs(eta,s.Rational(1,2))),('N_B_eta>0',s.diff(N_B,eta)),
 ('below A_0_eta implication',A_eta_num-s.Rational(7,2)*N_B.subs(eta,1-r)),('below S_eta>0',num(r-s_2)),
 ('below T_eta<0',num(-s.diff(T,eta)*eta**2)),('below T(.5)<r',num(r-T.subs(eta,s.Rational(1,2)))),
 ('N_2(etaB)<0 core',-G_B)]:check(label,expr,{r:rx,rho:px},[t_1,t_2])
check('below S(.5)>=r',num(S.subs(eta,s.Rational(1,2))-r),{r:rx,rho:px},[t_1,t_2],weak=True)

@dataclass(frozen=True)
class J:
    low:F;high:F
    @staticmethod
    def c(z):return z if isinstance(z,J) else J(F(z),F(z))
    def __add__(self,z):
        z=J.c(z);return J(self.low+z.low,self.high+z.high)
    __radd__=__add__
    def __neg__(self):return J(-self.high,-self.low)
    def __sub__(self,z):return self+-J.c(z)
    def __rsub__(self,z):return J.c(z)+-self
    def __mul__(self,z):
        z=J.c(z);rho=[self.low*z.low,self.low*z.high,self.high*z.low,self.high*z.high];return J(min(rho),max(rho))
    __rmul__=__mul__
root=json.loads((OUT/'local-root-certificate.json').read_text())
rad=F(root['radius']);rs,ps=[J(F(c)-rad,F(c)+rad) for c in root['center'][:2]]
def peval(expr,syms,box):
    poly=s.Poly(expr,*syms);powers=[]
    for value,degree in zip(box,poly.degree_list()):
        row=[J.c(1)]
        for _ in range(degree):row.append(row[-1]*value)
        powers.append(row)
    total=J.c(0)
    for mon,c in poly.terms():
        term=J.c(F(c))
        for row,j in zip(powers,mon):term=term*row[j]
        total=total+term
    return total
left,right=F(55278962,10**8),F(55278965,10**8)
jl=peval(N_1,[r,rho,eta],[rs,J.c(F(1,2)),J.c(left)])
jr=peval(N_1,[r,rho,eta],[rs,J.c(F(1,2)),J.c(right)])
je=peval(s.diff(N_1,eta),[r,rho,eta],[rs,J.c(F(1,2)),J(left,right)])
n2base=peval(N_2,[r,rho,eta],[rs,J.c(F(1,2)),J(left,right)])
assert jl.high<0 and jr.low>0 and je.low>0 and n2base.high<0
# Existing isolated root only: independently check the regularity used to
# infer R_resultant_rho=0 from H_root=0; no new root search or isolation is performed.
es=J(F(root['center'][2])-rad,F(root['center'][2])+rad)
curve_jac=s.diff(N_1,r)*s.diff(N_2,eta)-s.diff(N_1,eta)*s.diff(N_2,r)
regular=peval(curve_jac,[r,rho,eta],[rs,ps,es])
assert regular.low>0 or regular.high<0
coef=[]
for c in s.Poly(R_resultant,rho).all_coeffs():
    value=J.c(0)
    for a in s.Poly(c,r).all_coeffs():value=value*rs+F(a)
    coef.append(value)
for _ in range(2):
    out=[coef[0]]
    for c in coef[1:-1]:out.append(c+ps*out[-1])
    coef=out
# Interval Horner composition in rho=rs+(1-2rs)*x.
mapped=[J.c(0)];step=1-2*rs
for c in coef:
    nxt=[J.c(0)]*(len(mapped)+1)
    for i,a in enumerate(mapped):nxt[i]=nxt[i]+a*rs;nxt[i+1]=nxt[i+1]+a*step
    nxt[0]=nxt[0]+c;mapped=nxt
while len(mapped)>19:mapped.pop()
bern=[sum((mapped[j]*falling_ratio(k,j,18) for j in range(k+1)),J.c(0)) for k in range(19)]
assert all(c.high<0 for c in bern)
report={'independent_replay':True,'imports_producer':False,
        'method':'axis-by-axis Fraction power-to-Bernstein conversion; independent interval Horner deflation',
        'scope':'Full-rho kappa_1=0 dominant row face lemma chain; existing isolated root used',
        'positive_check_groups':len(checks),'positive_coefficient_count':sum(z['count'] for z in checks),
        'below_factor_identity_verified':True,'below_denominator_identity_verified':True,
        'left_identity_verified':True,
        'cap_expressions_reconstructed':True,'N_2_envelope_identity_verified':True,
        'H_root_tangent_identity_verified':True,'resultant_identity_verified':True,
        'base_lower_root_signs_verified':True,'stationary_curve_regularity_verified':True,
        'base_N_2_interval':[float(n2base.low),float(n2base.high)],
        'root_quotient_negative_coefficients':19,'all_checks_passed':True,
        'checks':checks}
(OUT/'full-rho-kink-independent-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:z for k,z in report.items() if k!='checks'},indent=2))
