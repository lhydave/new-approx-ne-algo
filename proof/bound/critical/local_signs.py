"""Exact local signs on the unique root box, using its rational kink identity."""
from pathlib import Path
from fractions import Fraction as F
import json,sys
import sympy as sp
OUT=Path(__file__).resolve().parent;BASE=OUT.parent.parent
sys.path.insert(0,str(BASE/'baseline_audit'))
from six_lp_shared_certificate import I
cert=json.loads((OUT/'local-root-certificate.json').read_text())
radius=F(cert['radius']);box=[I(F(c)-radius,F(c)+radius) for c in cert['center']]
r,rho,eta=box;sigma=1-rho;xi=sigma*(1-r)/(sigma+r)
lambda_0=r/rho;mu_0=r/sigma
vr=(1-eta-r)/eta;vc=(1-xi-r)/xi
s_1=1-lambda_0-(r-(1-eta)*(1-mu_0))/eta
s_2=1-mu_0-(r-(1-xi)*(1-lambda_0))/xi
h_2=1-mu_0+(rho/sigma)*vr;h_1=I.co(1)
S=1-(r-(1-eta)*s_2)/eta
T=h_2-(r-(1-xi)*s_1)/xi
alpha=(r-T)/(S-T);gamma=1-alpha
D_c=rho*(1-eta)+sigma*eta;D_r=sigma*(1-xi)+rho*xi
p_dagger=rho*(1-eta)/D_c;q_dagger=sigma*(1-xi)/D_r
J_C=rho*(r-eta*(1-lambda_0))/D_c;J_R=sigma*(r-xi*(1-mu_0))/D_r
base=(1-p_dagger)*(1-q_dagger)*r
knot=(base+q_dagger*(1-p_dagger*lambda_0)-J_R,base+p_dagger*(1-q_dagger*mu_0)-J_C)
left=((1-p_dagger)*r,(1-p_dagger)*r+p_dagger-J_C)
right=(1-p_dagger*lambda_0,p_dagger*(1-mu_0)-J_C)
below=((1-q_dagger)*r+q_dagger-J_R,(1-q_dagger)*r)
above=(q_dagger*(1-lambda_0)-J_R,1-q_dagger*mu_0)

def imax(x,y):return I(max(x.lo,y.lo),max(x.hi,y.hi))
def imin(vals):return I(min(x.lo for x in vals),min(x.hi for x in vals))
def edge(endpoint_a,endpoint_b):
    ar,ac=endpoint_a;br,bc=endpoint_b
    ar,ac,br,bc=map(I.co,[ar,ac,br,bc])
    candidates=[imax(ar,ac),imax(br,bc)]
    den=br-ar-bc+ac
    if den.lo>0 or den.hi<0:
        v=(ac-ar)/den
        if v.lo>=0 and v.hi<=1:candidates.append((1-v)*ar+v*br)
        else:assert v.hi<0 or v.lo>1,'crossing chart not isolated'
    else:assert False,'parallel crossing not isolated'
    return imin(candidates)

c={
 'r_lower':r-F(3,10),'r_upper':F(1,3)-r,
 'rho_lower':rho-r,'rho_upper':F(1,2)-rho,
 'eta_lower':eta-r,'eta_upper':1-r-eta,
 'xi_lower':xi-r,'xi_upper':1-r-xi,
 'alpha_lower':alpha-r,'alpha_upper':1-r-alpha,
 'gamma_lower':gamma-r,'gamma_upper':1-r-gamma,
 'first_row_source':eta*(1-lambda_0)+(1-eta)*(1-mu_0)-r,
 'first_col_source':xi*(1-mu_0)+(1-xi)*(1-lambda_0)-r,
 'column_deficit_range':eta*sigma/rho-(r-eta*(1-lambda_0)),
 'row_deficit_range':xi*rho/sigma-(r-xi*(1-mu_0)),
 'BC_positive':1-h_2,'S_minus_T_positive':S-T,
 'knot_row_nonnegative':knot[0],'knot_col_nonnegative':knot[1],
 'left_row_below_r':r-left[0],'left_col_above_r':left[1]-r,
 'knot_row_above_r':knot[0]-r,'knot_col_below_r':r-knot[1],
}
for name,z in [('vr',vr),('vc',vc),('s_1',s_1),('s_2',s_2),('h_2',h_2),('S',S),('T',T)]:
    c[name+'_positive']=z;c[name+'_below_one']=1-z
fr=[(0,1),(1-lambda_0,1-mu_0),(h_1,s_2)]
fc=[(1,0),(1-lambda_0,1-mu_0),(s_1,h_2)]
for name,pts in [('first_row',fr),('first_col',fc)]:
    for i,endpoint_a in enumerate(pts):
        for j,endpoint_b in enumerate(pts[i+1:],i+1):c[f'{name}_edge_{i}_{j}']=edge(endpoint_a,endpoint_b)-r
for i,(R,C) in enumerate([(vr,1),(s_1,h_2),(S,T),(1,0)]):
    if i!=2:c[f'alpha_source_{i}']=alpha*R+(1-alpha)*C-r
for i,(R,C) in enumerate([(1,vc),(h_1,s_2),(S,T),(0,1)]):
    if i!=2:c[f'gamma_source_{i}']=(1-gamma)*R+gamma*C-r
for name,endpoint_a in [('right',right),('below',below),('above',above)]:
    c['knot_'+name]=edge(endpoint_a,knot)-r
c['same_fixed_row_outer_chord']=edge(left,right)-r
c['same_fixed_col_outer_chord']=edge(below,above)-r
u=(1-eta)/eta;E=1-alpha+alpha*u
# The two gamma-basis terms of E_2, before its common threshold subtraction.
E_2_gamma_0_term=(1-alpha)*S+alpha*u*h_1
E_2_gamma_1_term=(1-alpha)*T+alpha*(1-r/eta+u*s_2)-r
c['column_envelope']=(1-gamma)*E_2_gamma_0_term+gamma*E_2_gamma_1_term-r*E
v=(1-xi)/xi
KR=1-r/xi+v*s_1
rowNR=(1-gamma)*S+gamma*KR-r
rowNC=(1-gamma)*T+gamma*v*h_2
c['row_envelope_decreases_alpha']=rowNC-rowNR
c['row_envelope_increases_gamma']=alpha*(KR-S)+(1-alpha)*(v*h_2-T)-r*(v-1)
assert all(z.lo>0 for z in c.values())

# Rational implicit-curve curvature at its stationary point.
r,rho,eta=sp.symbols('r rho eta');syms=[r,rho,eta]
meta=json.loads((OUT/'kink-polynomials.json').read_text())
n1,n2=[sp.sympify(meta[k]) for k in ['N_1','N_2']]
def peval(poly):
    terms=sp.Poly(poly,*syms).terms()
    deg=[max(m[k] for m,p_dagger in terms) for k in range(3)]
    powers=[]
    for val,d in zip(box,deg):
        ps=[I.co(1)]
        for _ in range(d):ps.append(ps[-1]*val)
        powers.append(ps)
    out=I.co(0)
    for mon,coef in terms:
        term=I.co(F(coef))
        for k,power in enumerate(mon):term=term*powers[k][power]
        out=out+term
    return out
ep=-peval(sp.diff(n1,rho))/peval(sp.diff(n1,eta))
quad=[]
for n in [n1,n2]:
    quad.append(peval(sp.diff(n,rho,2))+2*peval(sp.diff(n,rho,eta))*ep+peval(sp.diff(n,eta,2))*ep*ep)
den=peval(sp.diff(n1,r))*peval(sp.diff(n2,eta))-peval(sp.diff(n2,r))*peval(sp.diff(n1,eta))
rpp=-(quad[0]*peval(sp.diff(n2,eta))-quad[1]*peval(sp.diff(n1,eta)))/den
assert rpp.hi<F(-1)
report={'scope':'Exact signs inside the local unique-root box only; no global domination claim',
        'active_identities':['raw_h_1=1 / kappa_1=0','alpha_source_2=gamma_source_2=0','knot_left=0 via N_1=0','row_envelope=0 via N_2=0'],
        'all_nonactive_strictly_positive':True,'nonactive_count':len(c),
        'smallest_nonactive_residual_lower':float(min(z.lo for z in c.values())),
        'numeric_interval_signs':{key:[float(z.lo),float(z.hi)] for key,z in c.items()},
        'r_curve_second_derivative_interval':[float(rpp.lo),float(rpp.hi)],
        'r_curve_second_derivative_exact_upper':'-1'}
(OUT/'local-signs.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({key:val for key,val in report.items() if key!='numeric_interval_signs'},indent=2))
