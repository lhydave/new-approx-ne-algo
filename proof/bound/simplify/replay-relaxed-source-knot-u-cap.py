"""Replay the initial source-to-N_C orientation, source+L odds cap, and mirror signs."""
import sympy as s
from pathlib import Path
from fractions import Fraction
from math import comb
from itertools import product
import json

folder=Path(__file__).parent
saved=json.loads((folder.parent/'kink/u_cap_signs.json').read_text())
r,rho,u,v,a,b=s.symbols('r rho u v a b');sigma=1-rho;c=1-2*r
chi_1=1-r-r/rho;chi_2=1-r-r/sigma
N_R=sigma*sigma*chi_2+rho*(sigma-2*r)*u+sigma*sigma*c*v+rho*sigma*(chi_1-r)*u*v
N_C=rho*rho*chi_1+sigma*(rho-2*r)*v+rho*rho*c*u+rho*sigma*(chi_2-r)*u*v
L=(chi_1+c*u)*N_R+r*u*N_C
B_mix=(chi_2+c*v)*N_C+r*v*N_R
s_1=chi_1+chi_2*u;s_2=chi_2+chi_1*v
h_1_raw=1-r*(1+sigma)/rho+sigma*(1-r)*v/rho
h_2_raw=1-r*(1+rho)/sigma+rho*(1-r)*u/sigma
S=h_1_raw-r+u*(s_2-r);T=h_2_raw-r+v*(s_1-r)
uc=s.Rational(99,100)
sv=s.diff(S.subs(u,uc),v)
vs=-((S.subs(u,uc)-r).subs(v,0))/sv
J=rho*r-rho-199*r+100
assert s.factor(J-100*rho*sv)==0
assert s.factor(sigma*(1-r)/rho+chi_1-c/rho)==0
exprs={
 'K00':L.subs({u:uc,v:vs}),
 'K10':s.diff(L,u).subs({u:uc,v:vs}),
 'K01':s.diff(L,v).subs(u,uc),
 'K20':s.diff(L,u,2).subs(v,vs)/2,
 'K11':s.diff(L,u,v).subs(u,uc),
 'K21':s.diff(L,u,2,v)/2}
du,dv=s.symbols('du dv')
assert s.factor(L.subs({u:uc+du,v:vs+dv})-(exprs['K00']+exprs['K10']*du+exprs['K01']*dv+exprs['K20']*du**2+exprs['K11']*du*dv+exprs['K21']*du**2*dv))==0
# Initial orientation: under u>=uc and S>=r, the original source
# bounds force N_C<0. These seven checks precede use of L>=0.
S0=h_1_raw-r
Su=s_2-r
C0=N_C.subs(u,0)
Cu=s.diff(N_C,u)
assert s.factor(N_C-C0-u*Cu)==0
assert s.factor(s.diff(N_C,v)-sigma*(rho-2*r)-rho*sigma*(chi_2-r)*u)==0
assert s.factor(S.subs({u:uc,v:vs})-r)==0
CSraw=s.cancel((-Su)*C0+Cu*(S0-r))
CSclip=s.cancel((-Su)*C0+Cu*c)
assert s.cancel(N_C.subs(u,(S0-r)/(-Su))-CSraw/(-Su))==0
assert s.cancel(N_C.subs(u,c/(-Su))-CSclip/(-Su))==0
v0=r*(1+sigma)/(sigma*(1-r))
uh=(1-r)/r
def quadratic_coefficients(expr,left,right):
    assert s.Poly(expr,v).degree()<=2
    coeffs=[s.cancel(expr.subs(v,left)),
            s.cancel(expr.subs(v,left)+(right-left)*s.diff(expr,v).subs(v,left)/2),
            s.cancel(expr.subs(v,right))]
    t=s.symbols('t')
    reconstructed=coeffs[0]*(1-t)**2+2*coeffs[1]*t*(1-t)+coeffs[2]*t*t
    assert s.cancel(expr.subs(v,left+(right-left)*t)-reconstructed)==0
    return coeffs
raw_coeffs=quadratic_coefficients(CSraw,vs,v0)
clip_coeffs=quadratic_coefficients(CSclip,v0,uh)
initial_exprs=dict(Cuc0=s.cancel(N_C.subs({u:uc,v:vs})),
                  CSraw0=raw_coeffs[0],CSrawMid=raw_coeffs[1],CSrawEnd=raw_coeffs[2],
                  CSclip0=clip_coeffs[0],CSclipMid=clip_coeffs[1],CSclipEnd=clip_coeffs[2])
initial_denoms=dict(Cuc0=100*sigma*J,CSraw0=rho*sigma*J**2,
                   CSrawMid=2*rho*sigma*(1-r)*J,CSrawEnd=rho*sigma*(1-r)**2,
                   CSclip0=rho*sigma*(1-r)**2,CSclipMid=2*rho*r*sigma*(1-r),CSclipEnd=rho*r*r*sigma)
for name,denom in initial_denoms.items():
    recorded=s.sympify(saved[name]['positive_denominator'],locals={'r':r,'rho':rho})
    assert s.factor(recorded-denom)==0
exprs.update(initial_exprs)
rr=s.Rational(619,2000)+(s.Rational(1,3)-s.Rational(619,2000))*a
pp=rr+(s.Rational(1,2)-rr)*b
total=0
initial_total=0
for name,expr in exprs.items():
    cert=saved[name]
    denom=s.sympify(cert['positive_denominator'],locals={'r':r,'rho':rho})
    numerator=s.factor(expr*denom)
    assert s.denom(numerator)==1
    poly=s.Poly(s.expand(numerator.subs({r:rr,rho:pp})),a,b)
    powers={ks:Fraction(cv) for ks,cv in poly.terms()}
    da,db=cert['degree'];co=[]
    for i in range(da+1):
        row=[]
        for j in range(db+1):
            val=Fraction(0)
            for (h,k),cv in powers.items():
                if h<=i and k<=j:val+=cv*Fraction(comb(i,h),comb(da,h))*Fraction(comb(j,k),comb(db,k))
            assert val==Fraction(cert['coefficients'][i][j]) and val<0
            row.append(val);total+=1
            if name in initial_exprs:initial_total+=1
        co.extend(row)
    assert max(co)==Fraction(cert['max'])

mirror={rho:sigma,u:v,v:u}
for original,swapped in [(s_1,s_2),(s_2,s_1),(h_1_raw,h_2_raw),(h_2_raw,h_1_raw),(S,T),(T,S),(N_R,N_C),(N_C,N_R),(L,B_mix)]:
    assert s.factor(original.subs(mirror,simultaneous=True)-swapped)==0

# The both-high calculation uses no lower cap sign.  Its raw signs remain valid
# even when s_2<0: Y>=0, 0<=u<=1 suffice for T_v<0 and S_v>0.
assert s.factor(s.diff(T,v)-(s_1-r))==0
assert s.factor(s.diff(S,u)-(s_2-r))==0
assert s.factor(s.diff(S,v)-(sigma*(1-r)/rho+u*chi_1))==0
assert s.factor(chi_1+chi_2-(2*(1-r)-r/(rho*sigma)))==0
# rho*sigma<=1/4 gives s_1<=chi_1+chi_2<=2-6r, hence T_v<=2-7r<0.

both=json.loads((folder/'both-high-exclusion.json').read_text())
for chart in both['charts']:
    rho,r,z,w,t=s.symbols('rho r z w t')
    names={'rho':rho,'r':r,'u':u,'v':v,'z':z,'w':w,'t':t}
    pol=s.sympify(chart['polynomial'],locals=names)
    variables=[names[k] for k in chart['variables']]
    mapping={names[k]:s.sympify(value,locals=names) for k,value in chart['mapping'].items()}
    pp=s.Poly(s.expand(pol.subs(mapping)),*variables)
    ns=tuple(chart['degree']);powers={ks:Fraction(cv) for ks,cv in pp.terms()}
    for key,value in chart['coefficients'].items():
        ids=tuple(map(int,key.split(',')));val=Fraction(0)
        for ks,cv in powers.items():
            if not all(k<=i for k,i in zip(ks,ids)):continue
            factor=Fraction(1)
            for i,k,n in zip(ids,ks,ns):factor*=Fraction(comb(i,k),comb(n,k))
            val+=cv*factor
        assert val==Fraction(value) and val>0

report={'source_plus_L_only':True,'lower_cap_bounds_used':False,
 'six_shift_coefficients_strictly_negative':True,'negative_coefficients_recomputed':total,
 'seven_initial_N_C_checks_strictly_negative':True,
 'initial_N_C_negative_coefficients_recomputed':initial_total,
 'initial_N_C_source_substitution_and_quadratic_identities_verified':True,
 'denominator_common_factor':'J=100*rho*Sraw_v(99/100)>0',
 'all_mirror_identities_verified':True,'mirror_left_knot_equals_original_B_mix':True,
 'both_high_75_positive_coefficients_recomputed':True,
 'application':'Before shared-mass equalization: for rho<=1/2, S>=r and u>=99/100 force N_C<0 from original source bounds; the original left segment then gives L>=0, contradicting the six negative shift coefficients. The cap-only both-high lemma excludes simultaneous S,T>=r. Mirror identities are checked separately.',
 'result':'passed'}
(folder/'relaxed-source-knot-u-cap-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
