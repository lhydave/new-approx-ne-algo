"""Exact rational interval exclusion for the six-full-side-LP shared cap hulls.

All variables are free. The proof interval lower endpoint is selected after the
free diagnostic; it is not an algorithm parameter. Incomplete runs are NOT proofs.
"""
from fractions import Fraction as F
from dataclasses import dataclass
from pathlib import Path
import json,time,sys

@dataclass(frozen=True)
class I:
    lo:F
    hi:F
    def __post_init__(self):assert self.lo<=self.hi
    @staticmethod
    def co(x):return x if isinstance(x,I) else I(F(x),F(x))
    def __add__(self,b):
        b=I.co(b);return I(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,b):return self+-I.co(b)
    def __rsub__(self,b):return I.co(b)+-self
    def __mul__(self,b):
        b=I.co(b);p=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi];return I(min(p),max(p))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=I.co(b);assert b.lo>0 or b.hi<0
        return self*I(1/b.hi,1/b.lo)
    def __rtruediv__(self,b):return I.co(b)/self

def mx(x,c=0):return I(max(x.lo,F(c)),max(x.hi,F(c)))
def mn(x,c=1):return I(min(x.lo,F(c)),min(x.hi,F(c)))
def hull(vals):return I(min(vals),max(vals))
def repair_upper(S,T):
    if S.hi<0 or T.hi<0:return None
    a=min(F(1),S.hi);b=min(F(1),T.hi)
    assert 0<=a<=1 and 0<=b<=1
    return max(a,b)/(1+abs(a-b))

def edge_upper(a,b):
    # Coordinatewise upper caps are valid for the whole segment, so this exact
    # scalar minimax is an upper bound for the full fixed-side LP.
    ar,ac=[min(F(1),x.hi) for x in a]
    br,bc=[min(F(1),x.hi) for x in b]
    if min(ar,ac,br,bc)<0:return None
    out=min(max(ar,ac),max(br,bc))
    den=br-ar-bc+ac
    if den:
        t=(ac-ar)/den
        if 0<=t<=1:out=min(out,(1-t)*ar+t*br)
    return out

def caps_hull_upper(points):
    vals=[edge_upper(a,b) for i,a in enumerate(points) for b in points[i+1:]]
    return min(v for v in vals if v is not None) if any(v is not None for v in vals) else None

def evaluate(box,low):
    # Unit cube parameterization of the complete required domain.
    U,V,W,X=[I(a,b) for a,b in box]
    r=I(F(low)+(F(1,2)-F(low))*U.lo,F(low)+(F(1,2)-F(low))*U.hi)
    # rho is monotone in r and its own cube coordinate.
    rho=I(r.lo+(F(1,2)-r.lo)*V.lo,r.hi+(F(1,2)-r.hi)*V.hi)
    sigma=1-rho
    eta=hull([r+(1-2*r)*v for r in [r.lo,r.hi] for v in [W.lo,W.hi]])
    xi=hull([r+(1-2*r)*v for r in [r.lo,r.hi] for v in [X.lo,X.hi]])
    lambda_0=r/rho;mu_0=r/sigma
    D_c=rho+(sigma-rho)*eta;D_r=sigma+(rho-sigma)*xi
    # Sum representation ensures strictly positive denominators even on wide boxes.
    D_c=rho*(1-eta)+sigma*eta;D_r=sigma*(1-xi)+rho*xi
    p_dagger=rho*(1-eta)/D_c;q_dagger=sigma*(1-xi)/D_r
    lc=r-eta*(1-r/rho);lr=r-xi*(1-r/sigma)
    J_C=rho*lc/D_c;J_R=sigma*lr/D_r
    base=(1-p_dagger)*(1-q_dagger)*r
    mixr=base+q_dagger*(1-p_dagger*lambda_0)-J_R
    mixc=base+p_dagger*(1-q_dagger*mu_0)-J_C
    vr=mn((1-eta-r)/eta);vc=mn((1-xi-r)/xi)
    s_1=1-lambda_0-(r-(1-eta)*(1-mu_0))/eta
    s_2=1-mu_0-(r-(1-xi)*(1-lambda_0))/xi
    lower_rwj=mx(r/rho-(sigma/rho)*vc)
    lower_crz=mx(r/sigma-(rho/sigma)*vr)
    pairr=1-lower_rwj-(r-(1-eta)*mn(s_2))/eta
    pairc=1-lower_crz-(r-(1-xi)*mn(s_1))/xi
    constraints={
        'source_wz_w':eta*(1-lambda_0)+(1-eta)*(1-mu_0)-r,
        'source_wz_z':xi*(1-mu_0)+(1-xi)*(1-lambda_0)-r,
        'column_deficit_range':eta*sigma/rho-lc,
        'row_deficit_range':xi*rho/sigma-lr,
        'new_row_z_cap_nonnegative':s_1,
        'new_column_w_cap_nonnegative':s_2,
        'cross_row_cap_nonnegative':pairr,
        'cross_column_cap_nonnegative':pairc,
        'balanced_row_cap_nonnegative':mixr,
        'balanced_column_cap_nonnegative':mixc,
    }
    for name,val in constraints.items():
        if val.hi<0:return name,val.hi
    pr=repair_upper(pairr,pairc);mr=repair_upper(mixr,mixc)
    if pr is not None and pr<r.lo:return 'dual_pair_repair',pr-r.lo
    if mr is not None and mr<r.lo:return 'balanced_pair_repair',mr-r.lo
    h_1=1-lambda_0+(sigma/rho)*vc;h_2=1-mu_0+(rho/sigma)*vr
    Z=I.co(0);O=I.co(1)
    known_hulls={
      'first_row_shared_hull':[(Z,O),(1-lambda_0,1-mu_0),(h_1,s_2)],
      'first_column_shared_hull':[(O,Z),(1-lambda_0,1-mu_0),(s_1,h_2)],
      'second_row_shared_hull':[(vr,O),(s_1,h_2),(pairr,pairc),(O,Z)],
      'second_column_shared_hull':[(O,vc),(h_1,s_2),(pairr,pairc),(Z,O)],
    }
    for name,points in known_hulls.items():
        up=caps_hull_upper(points)
        if up is not None and up<r.lo:return name,up-r.lo
    return None

def encode_box(box):return [[str(a),str(b)] for a,b in box]
def run(low,node_limit):
    start=time.monotonic();stack=[(tuple((F(0),F(1)) for _ in range(4)),0,'')]
    leaves=[];nodes=0;hist={};maxdepth=0
    while stack and nodes<node_limit:
        box,depth,path=stack.pop();nodes+=1;maxdepth=max(maxdepth,depth)
        result=evaluate(box,low)
        if result:
            name,upper=result;hist[name]=hist.get(name,0)+1
            leaves.append({'path':path,'reason':name,'negative_upper':str(upper)})
        else:
            # Longest normalized side; tie broken cyclically through cube axes.
            widths=[b-a for a,b in box];axis=max(range(4),key=lambda k:(widths[k],-(k-depth)%4))
            a,b=box[axis];mid=(a+b)/2
            lower=list(box);upper=list(box);lower[axis]=(a,mid);upper[axis]=(mid,b)
            stack.append((tuple(upper),depth+1,path+str(axis)+'1'))
            stack.append((tuple(lower),depth+1,path+str(axis)+'0'))
        if nodes%5000==0:print(json.dumps({'nodes':nodes,'pending':len(stack),'leaves':len(leaves),'depth':maxdepth,'seconds':round(time.monotonic()-start,1),'hist':hist}),flush=True)
    report={'complete':not stack,'lower_bound':str(low),'upper_bound':'1/2','nodes':nodes,'leaves':len(leaves),'pending':len(stack),'max_depth':maxdepth,'seconds':time.monotonic()-start,'reasons':hist,'scope':'Exact rational necessary six-full-LP shared cap-hull model; no numerical samples used in exclusions'}
    out=Path(__file__).with_name('six-lp-shared-interval-'+str(low).replace('/','-')+'.json')
    out.write_text(json.dumps({'report':report,'leaves':leaves,'pending_boxes':[encode_box(b) for b,d,p in stack]},indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    low=F(sys.argv[1]) if len(sys.argv)>1 else F(8,25)
    limit=int(sys.argv[2]) if len(sys.argv)>2 else 100000
    run(low,limit)
