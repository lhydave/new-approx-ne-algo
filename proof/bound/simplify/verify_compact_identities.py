"""Small exact symbolic replay for the analytic simplification; no sampling."""
import json
from fractions import Fraction as F
from pathlib import Path
import sympy as s

rho, r, u, v, alpha, gamma, kappa_1, kappa_2 = s.symbols('rho r u v alpha gamma kappa_1 kappa_2')
sigma=1-rho
eta=1/(1+u); xi=1/(1+v)
chi_1=1-r-r/rho; chi_2=1-r-r/sigma
D_c_odds=rho*u+sigma; D_r_odds=sigma*v+rho
D_c=rho*(1-eta)+sigma*eta; D_r=sigma*(1-xi)+rho*xi
s_1=(chi_1+r)-r+u*((chi_2+r)-r); s_2=(chi_2+r)-r+v*((chi_1+r)-r)
h_1=1-kappa_1; h_2=1-kappa_2
S=h_1-r+u*(s_2-r); T=h_2-r+v*(s_1-r)
p_dagger=rho*u/D_c_odds; q_dagger=sigma*v/D_r_odds
J_C=rho*(r*(1+u)-(chi_1+r))/D_c_odds; J_R=sigma*(r*(1+v)-(chi_2+r))/D_r_odds
R=(1-p_dagger)*(1-q_dagger)*r+q_dagger*(1-p_dagger*r/rho)-J_R
C=(1-p_dagger)*(1-q_dagger)*r+p_dagger*(1-q_dagger*r/sigma)-J_C

checks={}
def check(name, lhs, rhs):
    assert s.factor(lhs-rhs)==0, name
    checks[name]=True

check('source_row', eta*(chi_1+r)+(1-eta)*(chi_2+r)-r, eta*s_1)
check('source_column', xi*(chi_2+r)+(1-xi)*(chi_1+r)-r, xi*s_2)
check('source_Dc', eta*s_1, 1-r-r*D_c/(rho*sigma))
check('source_Dr', xi*s_2, 1-r-r*D_r/(rho*sigma))
check('source_S', eta*h_1+(1-eta)*s_2-r, eta*S)
check('source_T', xi*h_2+(1-xi)*s_1-r, xi*T)
check('deficit_eta', eta*sigma/rho-(r-eta*(chi_1+r)), eta*(1-r)/rho-r)
check('deficit_xi', xi*rho/sigma-(r-xi*(chi_2+r)), xi*(1-r)/sigma-r)
check('BR_raw', r/rho-(sigma/rho)*((1-r)*v-r), (r*(1+sigma)-sigma*(1-r)*v)/rho)
check('BC_raw', r/sigma-(rho/sigma)*((1-r)*u-r), (r*(1+rho)-rho*(1-r)*u)/sigma)
check('knot_R', R, sigma*(((chi_2+r)-r+(1-r)*v)*D_c_odds+rho*r-r*u*v)/(D_c_odds*D_r_odds))
check('knot_C', C, rho*(((chi_1+r)-r+(1-r)*u)*D_r_odds+sigma*r-r*u*v)/(D_c_odds*D_r_odds))
check('left_R', (1-p_dagger)*r, sigma*r/D_c_odds)
check('left_C', (1-p_dagger)*r+p_dagger-J_C, rho*((1-r)*u+1-2*r)/D_c_odds)
check('right_R', 1-p_dagger*r/rho, 1-r*u/D_c_odds)
check('right_C', p_dagger*(chi_2+r)-J_C, rho*s_1/D_c_odds)
check('below_R', (1-q_dagger)*r+q_dagger-J_R, sigma*((1-r)*v+1-2*r)/D_r_odds)
check('below_C', (1-q_dagger)*r, rho*r/D_r_odds)
check('above_R', q_dagger*(chi_1+r)-J_R, sigma*s_2/D_r_odds)
check('above_C', 1-q_dagger*r/sigma, 1-r*v/D_r_odds)

H=alpha*S+(1-alpha)*T; K=(1-gamma)*S+gamma*T
U=alpha*s_1+(1-alpha)*h_2; V=(1-gamma)*h_1+gamma*s_2
old_row=alpha*((1-gamma)*S+gamma*(1-r/xi+v*s_1)-r)+(1-alpha)*((1-gamma)*T+gamma*v*h_2)-r*(1-gamma+gamma*v)
old_col=(1-gamma)*((1-alpha)*S+alpha*u*h_1)+gamma*((1-alpha)*T+alpha*(1-r/eta+u*s_2)-r)-r*(1-alpha+alpha*u)
check('row_surplus', old_row, (1-gamma)*(H-r)+gamma*v*(U-r)+alpha*gamma*(1-r-r*v)-alpha*r)
check('column_surplus', old_col, (1-alpha)*(K-r)+alpha*u*(V-r)+alpha*gamma*(1-r-r*u)-gamma*r)
df=rho+r*(sigma-rho)
check('above_quadratic', (r-(sigma*(1-r)/df-r/rho))*rho*df,
      (1-2*r*r)*rho*rho-(1-r)**2*rho+r*r)

assert F(21,25)>F(361,441)
assert F(11,35)<F(1,3)
assert F(19,30)-F(2,3)==-F(1,30)
assert F(37,175)<F(3,10)
assert F(2401,10000)<F(7,25)
checks['exact_bound_constants']=True
report={'scope':'Exact symbolic identities and rational constants only; global inequality proofs in companion derivation', 'checks':checks, 'all_passed':True}
out=Path(__file__).with_name('compact-identities-verification.json')
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'all_passed':True,'checks':len(checks),'output':str(out)}))
