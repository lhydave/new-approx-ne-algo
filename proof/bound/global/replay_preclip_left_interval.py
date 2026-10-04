"""Independent high-rho L=0 binary-cover replay; imports no producer.

The interval implementation comes from our independent B_mix checker. The L
model, guarded contractors and cubic Bernstein identity are reconstructed here.
No tree search, floating point optimizer or production evaluate() is imported.
"""
from fractions import Fraction as ExactFraction
from pathlib import Path
from collections import Counter
import argparse,json,time,hashlib
from replay_preclip_below_interval import Range

LOW=ExactFraction(30953996,10**8);HIGH=ExactFraction(1,3)

def independent_contradiction(box):
    ar,ap,au=(Range(lo,hi)for lo,hi in box)
    r=LOW+(HIGH-LOW)*ar
    rho=ExactFraction(1,2)+(ExactFraction(1,2)-r)*ap;sigma=1-rho
    assert rho.lower>0 and sigma.lower>0 and (1-r).lower>0
    k=r/(1-r);u=k+(1-k)*au;c=1-2*r
    end=k+r/(sigma*(1-r))
    # sigma>=r and r<=1/3 imply r(1+sigma)<=2q(1-r), hence end<=2.
    v=Range(k.lower,min(end.upper,ExactFraction(2)))
    assert v.lower>0
    chi_1=1-r-r/rho;chi_2=1-r-r/sigma
    row0=sigma*sigma*chi_2+rho*(sigma-2*r)*u
    rowv=sigma*sigma*c+rho*sigma*(chi_1-r)*u
    col0=rho*rho*chi_1+rho*rho*c*u
    colv=sigma*(rho-2*r)+rho*sigma*(chi_2-r)*u
    numerator=sigma*chi_2*(sigma+rho*u)*(chi_1+c*u)
    denominator=-((chi_1+c*u)*rowv+r*u*colv)
    # L=numerator-denominator*v. No global sign on denominator is used.
    for repeat in range(3):
        implied=numerator/v
        denominator=denominator.intersect(implied.lower,implied.upper)
        if denominator is None:return 'affine-root-empty'
        if denominator.lower>0 or denominator.upper<0:
            implied=numerator/denominator
            v=v.intersect(implied.lower,implied.upper)
            if v is None:return 'affine-root-empty'
        residual=numerator-denominator*v
        if residual.lower>0 or residual.upper<0:return 'affine-root-empty'

    h_2=1-r*(1+rho)/sigma+rho*(1-r)*u/sigma
    s_1=chi_1+chi_2*u
    gap_r=r-s_1;gap_h=h_2-r
    if gap_r.upper<=0:return 's_1>=r'
    if gap_h.upper<=0:return 'h_2<=r'
    gap_r=gap_r.intersect(lower=0);gap_h=gap_h.intersect(lower=0)
    intercept_m=1-r*(1+sigma)/rho-2*r+(chi_2-r)*u
    slope_m=sigma*(1-r)/rho+chi_1*u
    intercept_n=2*r-h_2
    linear_budgets=[
        (intercept_m,slope_m),
        (intercept_n,gap_r),
        (chi_2,chi_1),
        (gap_h,-gap_r),
        ((1-r)*intercept_m-r*intercept_n,(1-r)*slope_m-r*gap_r),
        (gap_h*intercept_m-gap_r*intercept_n,gap_h*slope_m-gap_r.square())]
    # For f(v)=a+b*v>=0, a<=a_hi,b<=b_hi,v>0 gives the necessary
    # relaxation a_hi+b_hi*v>=0. Divisions use b_hi strictly nonzero.
    for repeat in range(4):
        for a,b in linear_budgets:
            if b.upper<0:
                v=v.intersect(upper=-a.upper/b.upper)
            elif b.upper>0 and a.upper<0:
                v=v.intersect(lower=-a.upper/b.upper)
            elif b.upper==0 and a.upper<0:
                return 'source-contract-empty'
            if v is None:return 'source-contract-empty'
        implied=numerator/v
        denominator=denominator.intersect(implied.lower,implied.upper)
        if denominator is None:return 'affine-root-empty'
        if denominator.lower>0 or denominator.upper<0:
            implied=numerator/denominator
            v=v.intersect(implied.lower,implied.upper)
            if v is None:return 'affine-root-empty'

    s_1=chi_1+chi_2*u;s_2=chi_2+chi_1*v
    if s_2.upper<0:return 's_2<0'
    if s_2.lower>=r.upper:return 's_2>=r'
    s_2=s_2.intersect(lower=0,upper=r.upper)
    h_1=1-r*(1+sigma)/rho+sigma*(1-r)*v/rho
    if h_1.lower>1:return 'preclip-domain-empty'
    h_1=h_1.intersect(upper=1)
    h_2=1-r*(1+rho)/sigma+rho*(1-r)*u/sigma
    S=h_1-r+(s_2-r)*u;T=h_2-r+(s_1-r)*v
    m=S-r;n=r-T
    if m.upper<=0:return 'm<=0'
    if n.upper<=0:return 'n<=0'
    if T.upper<0:return 'T<0'
    m=m.intersect(lower=0);n=n.intersect(lower=0,upper=r.upper)
    B_mass=(1-r)*m-r*n
    if B_mass.upper<0:return 'upper-mass<0'
    gap_r=r-s_1;gap_h=h_2-r
    if gap_r.upper<=0:return 's_1>=r'
    if gap_h.upper<=0:return 'h_2<=r'
    gap_r=gap_r.intersect(lower=0);gap_h=gap_h.intersect(lower=0)
    G=m*gap_h-n*gap_r
    if G.upper<0:return 'G<0'
    homogeneous=v*gap_h*m.square()+n*m*(c-v*(2*r-s_1))-r*n.square()
    if homogeneous.upper<0:return 'Ehom<0'

    # Keep the cubic v dependence, using independently expanded coefficients.
    b=2*r-s_1
    coeff=[
        c*intercept_n*intercept_m-r*intercept_n.square(),
        gap_h*intercept_m.square()+c*(intercept_n*slope_m+gap_r*intercept_m)
            -b*intercept_n*intercept_m-2*r*intercept_n*gap_r,
        2*gap_h*intercept_m*slope_m+c*gap_r*slope_m
            -b*(intercept_n*slope_m+gap_r*intercept_m)-r*gap_r.square(),
        gap_h*slope_m.square()-b*gap_r*slope_m]
    lo=v.lower;span=v.upper-v.lower
    power=[coeff[0]+coeff[1]*lo+coeff[2]*lo*lo+coeff[3]*lo*lo*lo,
        span*(coeff[1]+2*coeff[2]*lo+3*coeff[3]*lo*lo),
        span*span*(coeff[2]+3*coeff[3]*lo),span*span*span*coeff[3]]
    bernstein=[power[0],power[0]+power[1]/3,
        power[0]+2*power[1]/3+power[2]/3,sum(power,Range(0))]
    if all(z.upper<0 for z in bernstein):return 'Ehom-Bernstein<0'
    if m.lower>0:
        ratio=(n/m).intersect(lower=0,upper=(1-r.lower)/r.lower)
        if ratio is None:return 'upper-mass<0'
        envelope=v*gap_h+ratio*(c-v*(2*r-s_1))-r*ratio.square()
        if envelope.upper<0:return 'E<0'
    return None

def symbolic_audit():
    import sympy as s
    v,A,B,C,R,M,N,W=s.symbols('v A B C R M N W')
    expression=v*A*(M+W*v)**2+(N+R*v)*(M+W*v)*(C-B*v)-s.Symbol('r')*(N+R*v)**2
    r=s.Symbol('r')
    independently_expanded=[C*N*M-r*N**2,
        A*M**2+C*(N*W+R*M)-B*N*M-2*r*N*R,
        2*A*M*W+C*R*W-B*(N*W+R*M)-r*R**2,A*W**2-B*R*W]
    assert s.expand(expression-sum(a*v**j for j,a in enumerate(independently_expanded)))==0
    r,rho,u=s.symbols('r rho u');sigma=1-rho;c=1-2*r;chi_1=1-r-r/rho;chi_2=1-r-r/sigma
    N_R=sigma*sigma*chi_2+rho*(sigma-2*r)*u+(sigma*sigma*c+rho*sigma*(chi_1-r)*u)*v
    N_C=rho*rho*chi_1+rho*rho*c*u+(sigma*(rho-2*r)+rho*sigma*(chi_2-r)*u)*v
    L=s.expand((chi_1+c*u)*N_R+r*u*N_C)
    assert s.factor(L.subs(v,0)-sigma*chi_2*(sigma+rho*u)*(chi_1+c*u))==0
    assert s.diff(L,v,2)==0
    return True

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('certificate',nargs='?',default=str(Path(__file__).resolve().parent.parent/'critical/preclip-left-interval-certificate.json'))
    args=ap.parse_args();path=Path(args.certificate)
    assert path.stat().st_size<=10*1024*1024
    cert=json.loads(path.read_text())
    assert cert['version']==1 and not cert['upper_root_condition']
    assert ExactFraction(cert['r_lower'])==LOW and ExactFraction(cert['r_upper'])==HIGH
    assert cert['root_domain']=='[0,1]^3'
    assert symbolic_audit()
    nodes=cert['nodes'];assert len(nodes)<=199999
    stack=[(0,((ExactFraction(0),ExactFraction(1)),)*3,0)];seen=set()
    independent=Counter();recorded=Counter();closed=opened=0;depth_max=0
    start=time.monotonic()
    while stack:
        index,box,depth=stack.pop()
        assert isinstance(index,int) and 0<=index<len(nodes) and index not in seen
        seen.add(index);node=nodes[index];depth_max=max(depth_max,depth)
        if 'axis' in node:
            assert set(node)=={'axis','left','right'} and node['axis']in(0,1,2)
            axis=node['axis'];assert node['left']!=node['right']
            lo,hi=box[axis];middle=(lo+hi)/2
            left=list(box);right=list(box);left[axis]=(lo,middle);right[axis]=(middle,hi)
            stack.extend([(node['right'],tuple(right),depth+1),(node['left'],tuple(left),depth+1)])
        elif node.get('open')is True:
            assert set(node)=={'open'};opened+=1
        else:
            assert set(node)=={'reason'} and isinstance(node['reason'],str)
            reason=independent_contradiction(box)
            assert reason is not None,('unverified leaf',index,node['reason'],box)
            assert reason==node['reason'],('reason mismatch',index,node['reason'],reason)
            independent[reason]+=1;recorded[node['reason']]+=1;closed+=1
        assert time.monotonic()-start<900,'replay time budget'
    assert len(seen)==len(nodes) and closed+opened==cert['leaf_count']<=100000
    assert closed==cert['closed_leaves'] and opened==cert['open_leaves']==0
    assert dict(recorded)==cert['reasons'] and depth_max==cert['maximum_depth']
    assert cert['all_leaves_closed']
    assert 5*LOW**4+2*LOW**2-4*LOW+1<0
    assert 20*HIGH**3+4*HIGH-4<0
    report={'success':True,'scope':'high-rho all physical preclip L=0 roots; no upper-root restriction',
        'conditional_boundary_theorem_only':True,'independent_fraction_intervals':True,
        'exact_binary_cover':True,'safe_denominators':True,'pole_intervals_preserved':True,
        'affine_contractors_audited':True,'cubic_Bernstein_identity_verified':True,
        'tree_nodes':len(nodes),'closed_leaves_verified':closed,'open_leaves':opened,
        'maximum_depth':depth_max,'all_recorded_reasons_recomputed':True,
        'independent_reasons':dict(independent),'seconds':time.monotonic()-start,
        'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    target=Path(__file__).with_name('preclip-left-independent-replay.json')
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
