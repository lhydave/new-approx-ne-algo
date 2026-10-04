"""Symbolic identities only; no optimizer and no global-bound assertion."""
import sympy as s

r,alpha,gamma,u,v,S,T,h_2,h_1=s.symbols('r alpha gamma u v S T h_2 h_1')
s_1=(T-h_2+r*(1+v))/v
s_2=(S-h_1+r*(1+u))/u
Fr=alpha*((1-gamma)*S+gamma*(1-r*(1+v)+v*s_1)-r)+(1-alpha)*((1-gamma)*T+gamma*v*h_2)-r*(1-gamma+gamma*v)
Fc=(1-gamma)*((1-alpha)*S+alpha*u*h_1)+gamma*((1-alpha)*T+alpha*(1-r*(1+u)+u*s_2)-r)-r*(1-alpha+alpha*u)
d=S-T; A_0=(r-T)/d
assert s.factor(s.diff(Fr,alpha)-((1-gamma)*d+gamma*(1+T-h_2*(1+v))-r))==0
assert s.factor(Fr.subs({alpha:A_0,gamma:0})+A_0*r)==0
assert s.factor(s.diff(Fr,alpha).subs(gamma,r)-((1-r)*S-(1-2*r)*T-r*h_2*(1+v)))==0
endpoint=((S-r)*(1-h_2*(1+v))+T*(T-r))/d
assert s.factor(s.diff(Fr,alpha).subs(gamma,1-A_0)-endpoint)==0
n=r-T;m=S-r
redr=n*m*(1+T-h_2)+v*m*m*h_2-r*d*(n+v*m)
redc=n*m*(1+S-h_1)+u*n*n*h_1-r*d*(m+u*n)
assert s.factor(d*d*Fr.subs({alpha:A_0,gamma:1-A_0})-redr)==0
assert s.factor(d*d*Fc.subs({alpha:A_0,gamma:1-A_0})-redc)==0

rho,sigma=s.symbols('rho sigma');D_c_odds=rho*u+sigma;D_r_odds=sigma*v+rho
p_dagger=rho*u/D_c_odds;q_dagger=sigma*v/D_r_odds
J_C=(rho*r*(1+u)-rho+r)/D_c_odds;J_R=(sigma*r*(1+v)-sigma+r)/D_r_odds
Rk=(1-p_dagger)*(1-q_dagger)*r+q_dagger*(1-p_dagger*r/rho)-J_R
Ck=(1-p_dagger)*(1-q_dagger)*r+p_dagger*(1-q_dagger*r/sigma)-J_C
fR=r*(1+sigma)-sigma*(1-r)*v;fC=r*(1+rho)-rho*(1-r)*u
Rk2=(D_c_odds*(sigma-fR)+rho*sigma*r-sigma*u*v*r)/(D_c_odds*D_r_odds)
Ck2=(D_r_odds*(rho-fC)+rho*sigma*r-rho*u*v*r)/(D_c_odds*D_r_odds)
assert s.factor((Rk-Rk2).subs(sigma,1-rho))==0
assert s.factor((Ck-Ck2).subs(sigma,1-rho))==0
assert s.factor(((1-p_dagger)*r+p_dagger-J_C-r-(rho-fC-rho*u*r)/D_c_odds).subs(sigma,1-rho))==0
print('All 10 reduced-envelope and knot identities verified symbolically.')
