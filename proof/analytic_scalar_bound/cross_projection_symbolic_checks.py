from pathlib import Path
import sympy as s,json
r,g,rho,eta,xi,p,q,u,v=s.symbols('r g rho eta xi p q u v');sigma=1-rho
D_c=rho+(sigma-rho)*eta;D_r=sigma+(rho-sigma)*xi;p_dagger=rho*(1-eta)/D_c;q_dagger=sigma*(1-xi)/D_r
J_C=(rho*(r-eta)+eta*g)/D_c;J_R_slope=(r-xi*(1-g/sigma))/(1-xi);B=(1-p_dagger)*g;U=1-p_dagger*g/rho-J_R_slope;D=U+p_dagger*g/sigma
F=D*(B-r)+(p_dagger-J_C)*(U-B)
checks={}
def ck(name,expr):
 z=s.factor(expr);assert z==0,(name,z);checks[name]=True
D_c_odds=rho*u+sigma
N=sigma*(rho*(2-r-rho*(1-r)-2*g)-r)*u*u*v+2*sigma*sigma*(rho*(1-r)-g)*u*u+rho*(rho*rho*r-rho*rho-2*rho*g+rho*r+rho-3*r+1)*u*v+sigma*(2*rho*(1-2*g)-r)*u+rho*rho*(1-2*r)*v+rho*(3*rho*r-rho-2*r+1)
N=s.cancel(u**2*v*N.subs({u:1/u,v:1/v},simultaneous=True))
ck('odds_Fg_numerator',s.diff(F,g).subs({eta:1/(1+u),xi:1/(1+v)})+N/(sigma*v*D_c_odds**2))
ck('odds_Ng_factor',s.diff(N,g)+2*(rho*sigma+sigma**2*v+rho**2*u+2*rho*sigma*u*v))
Geta=rho*sigma*(1-r)/D_c;Gxi=rho*sigma*(1-r)/D_r
ck('source_upper_chart_order',Geta-Gxi+rho*sigma*(1-r)*(sigma-rho)*(eta+xi-1)/(D_c*D_r))
C0=B+p_dagger-J_C
ck('source_first_col_g_derivative',s.diff(C0,g)+rho*eta/D_c)
ck('source_first_row_slope_lower_difference',(U-B)-(1-p_dagger)*(g/rho-g)-(1-g/rho-J_R_slope))
R1=(1-p_dagger)*(1-q_dagger)*g+q_dagger*(1-p_dagger*g/rho)-(sigma*(r-xi)+xi*g)/D_r
C1=(1-p_dagger)*(1-q_dagger)*g+p_dagger*(1-q_dagger*g/sigma)-J_C
ck('joint_row_g_derivative',s.diff(R1,g)-(xi/D_r*(rho*(1-p_dagger)-1)-p_dagger*q_dagger/rho))
ck('joint_col_g_derivative',s.diff(C1,g)-(eta/D_c*(sigma*(1-q_dagger)-1)-p_dagger*q_dagger/sigma))
Fendpoint=(R1-r)*(C0-r)-(r-B)*(r-C1)
ck('crossing_determinant_divided_by_knot',Fendpoint-q_dagger*F)
Rhigh=(1-p_dagger)*(1-q)*g+q*(1-p_dagger*g/rho)-(1-q)*(r-xi*(1-g/sigma))*sigma/(xi*rho)
Chigh=(1-p_dagger)*(1-q)*g+p_dagger*(1-q*g/sigma)-J_C
ck('outer_knot_row_g_derivative',s.diff(Rhigh,g)-((1-q)*(1-p_dagger-1/rho)-p_dagger*q/rho))
ck('outer_knot_col_g_derivative',s.diff(Chigh,g)-((1-p_dagger)*((1-q)-1/sigma)-p_dagger*q/sigma))
s_1=1-g/rho-(r-(1-eta)*(1-g/sigma))/eta
s_2=1-g/sigma-(r-(1-xi)*(1-g/rho))/xi
ck('source_s_1_tight_g_derivative',s.diff(s_1,g)+1/rho+(1-eta)/(eta*sigma))
ck('source_s_2_tight_g_derivative',s.diff(s_2,g)+1/sigma+(1-xi)/(xi*rho))
# Independently reconstruct each certificate polynomial from simultaneous maps.
num,den=s.fraction(s.factor(s.diff(F,g)));K=-num
zvars=s.symbols('z_r z_rho z_eta z_xi');zr,zp,ze,zx=zvars
rmap=s.Rational(3,10)+zr/30;emap=r+(1-2*r)*ze
for label,large in [('cross_derivative_certificate',False),('cross_derivative_large_mass_certificate',True)]:
 raw=json.loads(Path(__file__).with_name(label+'.json').read_text())
 for chart in raw['charts']:
  ge=chart['chart']=='sum_ge1'
  expected=D_r if ge==large else D_c
  upper=rho*sigma*(1-r)/expected
  kn,kd=s.fraction(s.factor(K.subs(g,upper)));kn=s.factor(kn/(kd/expected))
  xmap=1-eta+(eta-r)*zx if ge else r+(1-eta-r)*zx
  Pmap=s.Rational(1,2)+(s.Rational(1,2)-r)*zp if large else r+(s.Rational(1,2)-r)*zp
  # Recompose maps into expressions in only the four independent cube variables.
  final_e=emap.subs(r,rmap);final_P=Pmap.subs(r,rmap)
  final_x=xmap.subs(eta,final_e).subs(r,rmap)
  direct=s.expand(kn.subs({r:rmap,rho:final_P,eta:final_e,xi:final_x},simultaneous=True))
  stored=s.sympify(chart['polynomial'],locals={str(z):z for z in zvars})
  ck(label+'_'+chart['chart']+'_mapped_identity',direct-stored)
result={'scope':'Independent symbolic checks for all monotonicity reductions and producer-to-cube polynomial identities.','passed':len(checks),'checks':checks}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
