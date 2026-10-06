"""Independent replay of the bounded conditional B_mix-boundary certificate.

Imports no proof producer. Rational interval arithmetic is implemented here;
the midpoint binary cover, every contraction, and every contradiction leaf are
recomputed. A partial tree never receives a theorem success flag.
"""
from fractions import Fraction as ExactFraction
from pathlib import Path
from collections import Counter
import json,time,argparse

class Range:
    def __init__(self,lo,hi=None):
        self.lower=ExactFraction(lo);self.upper=ExactFraction(lo if hi is None else hi)
        assert self.lower<=self.upper
    @staticmethod
    def cast(value):return value if isinstance(value,Range) else Range(value)
    def __add__(self,value):
        other=self.cast(value)
        return Range(self.lower+other.lower,self.upper+other.upper)
    __radd__=__add__
    def __neg__(self):return Range(-self.upper,-self.lower)
    def __sub__(self,value):return self+-self.cast(value)
    def __rsub__(self,value):return self.cast(value)+-self
    def __mul__(self,value):
        other=self.cast(value)
        corners=(self.lower*other.lower,self.lower*other.upper,
                 self.upper*other.lower,self.upper*other.upper)
        return Range(min(corners),max(corners))
    __rmul__=__mul__
    def __truediv__(self,value):
        other=self.cast(value)
        assert other.upper<0 or other.lower>0,'division through zero'
        return self*Range(1/other.upper,1/other.lower)
    def __rtruediv__(self,value):return self.cast(value)/self
    def square(self):
        small=ExactFraction(0) if self.lower<=0<=self.upper else min(self.lower**2,self.upper**2)
        return Range(small,max(self.lower**2,self.upper**2))
    def intersect(self,lower=None,upper=None):
        lo=self.lower if lower is None else max(self.lower,ExactFraction(lower))
        hi=self.upper if upper is None else min(self.upper,ExactFraction(upper))
        return Range(lo,hi) if lo<=hi else None

rlo=ExactFraction(30953996,10**8);rhi=ExactFraction(1,3)

def independent_exclusion(box):
    # Reconstruct the coupled physical map, rather than reading range endpoints
    # or cap/envelope values from the certificate.
    a,b,w=(Range(lo,hi) for lo,hi in box)
    r=rlo+(rhi-rlo)*a
    rho=ExactFraction(1,2)+(ExactFraction(1,2)-r)*b;sigma=1-rho
    assert rho.lower>0 and sigma.lower>0 and (1-r).lower>0
    k=r/(1-r);v=k+r*w/(sigma*(1-r));c=1-2*r
    chi_1=1-r-r/rho;chi_2=1-r-r/sigma
    below=chi_2+c*v
    # At an actual opposite B_mix=0 boundary, below>0, chi_1>0, and B0>0.
    if below.upper<=0:return 'below-outer<=r'
    if chi_1.upper<=0:return 'chi_1<=0'
    xpos=chi_1.intersect(lower=0);bpos=below.intersect(lower=0)
    B0=rho*xpos*(rho+sigma*v)*bpos
    row_u=rho*(sigma-2*r)+rho*sigma*(chi_1-r)*v
    col_u=rho*rho*c+rho*sigma*(chi_2-r)*v
    den=-(below*col_u+r*v*row_u)
    if den.upper<=0:return 'root-den<=0'
    den=den.intersect(lower=0)
    u=Range(k.lower,1)
    # B_mix=B0-den*u=0. Since u>=k>0 and den>0, both intersections
    # below are necessary, including when the box straddles den=0.
    for repeat in range(2):
        implied_den=B0/u
        den=den.intersect(implied_den.lower,implied_den.upper)
        if den is None:return 'root-contract-empty'
        if den.upper<=0:return 'root-den<=0'
        implied_low=B0.lower/den.upper
        implied_high=B0.upper/den.lower if den.lower>0 else ExactFraction(1)
        u=u.intersect(implied_low,implied_high)
        if u is None:return 'root-contract-empty'

    s_1=chi_1+chi_2*u;s_2=chi_2+chi_1*v
    if s_2.upper<0:return 's_2<0'
    if s_2.lower>=r.upper:return 's_2>=r'
    # These intersections use s_2>=0 and the independently proved high-range
    # cap s_2<r. They do not assume s_1>=0.
    s_2=s_2.intersect(lower=0,upper=r.upper)
    row=sigma*sigma*chi_2+row_u*u+sigma*sigma*c*v
    col=rho*rho*chi_1+sigma*(rho-2*r)*v+col_u*u
    if row.lower>=0:return 'N_R>=0'
    if col.upper<0:return 'N_C<0'
    h_1=1-r*(1+sigma)/rho+sigma*(1-r)*v/rho
    h_2=1-r*(1+rho)/sigma+rho*(1-r)*u/sigma
    S=h_1-r+(s_2-r)*u;T=h_2-r+(s_1-r)*v
    m=S-r;n=r-T
    if m.upper<=0:return 'm<=0'
    if n.upper<=0:return 'n<=0'
    if T.upper<0:return 'T<0'
    m=m.intersect(lower=0)
    n=n.intersect(lower=0,upper=r.upper)  # T>=0 gives n<=r.
    Theta=(1-r)*m-r*n
    if Theta.upper<0:return 'upper-mass<0'
    R=r-s_1;A=h_2-r
    if R.upper<=0:return 's_1>=r'
    # G>=0, m,n,R>0 imply A>0, so this intersection is conditional but sound.
    if A.upper<=0:return 'h_2<=r'
    R=R.intersect(lower=0);A=A.intersect(lower=0)
    G=m*A-n*R
    if G.upper<0:return 'G<0'
    Ehom=v*A*m.square()+n*m*(c-v*(2*r-s_1))-r*n.square()
    if Ehom.upper<0:return 'Ehom<0'
    if m.lower>0:
        z=(n/m).intersect(lower=0,upper=(1-r.lower)/r.lower)
        if z is None:return 'upper-mass<0'
        E=v*A+z*(c-v*(2*r-s_1))-r*z.square()
        if E.upper<0:return 'E<0'
    return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('certificate',nargs='?',default=str(Path(__file__).resolve().parent.parent/'critical/preclip-below-interval-certificate.json'))
    args=ap.parse_args();path=Path(args.certificate)
    assert path.stat().st_size<=10*1024*1024
    cert=json.loads(path.read_text());assert cert['version']==1
    assert ExactFraction(cert['r_lower'])==rlo and ExactFraction(cert['r_upper'])==rhi
    assert cert['root_domain']=='[0,1]^3'
    nodes=cert['nodes'];assert len(nodes)<=9999
    stack=[(0,((ExactFraction(0),ExactFraction(1)),)*3,0)];seen=set()
    independent_reasons=Counter();recorded_reasons=Counter();closed=opened=0
    maxdepth=0;start=time.monotonic()
    while stack:
        index,box,depth=stack.pop();assert index not in seen
        assert isinstance(index,int) and 0<=index<len(nodes)
        seen.add(index);node=nodes[index];maxdepth=max(maxdepth,depth)
        if 'axis' in node:
            assert set(node)=={'axis','left','right'}
            axis=node['axis'];assert axis in (0,1,2)
            left,right=node['left'],node['right'];assert left!=right
            lo,hi=box[axis];mid=(lo+hi)/2
            lbox=list(box);rbox=list(box)
            lbox[axis]=(lo,mid);rbox[axis]=(mid,hi)
            stack.extend([(right,tuple(rbox),depth+1),(left,tuple(lbox),depth+1)])
        elif node.get('open') is True:
            assert set(node)=={'open'};opened+=1
        else:
            assert set(node)=={'reason'} and isinstance(node['reason'],str)
            reason=independent_exclusion(box)
            assert reason is not None,('unverified leaf',index,node['reason'],box)
            independent_reasons[reason]+=1;recorded_reasons[node['reason']]+=1;closed+=1
        assert time.monotonic()-start<600,'independent replay budget exceeded'
    assert len(seen)==len(nodes),'unreachable certificate nodes'
    assert closed+opened==cert['leaf_count']<=5000
    assert closed==cert['closed_leaves'] and opened==cert['open_leaves']
    assert dict(recorded_reasons)==cert['reasons']
    assert maxdepth==cert['maximum_depth']
    assert cert['all_leaves_closed']==(opened==0)
    # Check the exact endpoint of the external uniform-cap lemma. Its derivative
    # is negative to 1/3, hence s_1,s_2<r is a genuine dependency, not a new claim
    # extracted from an interval box.
    assert 5*rlo**4+2*rlo**2-4*rlo+1<0
    report={'success':opened==0,'conditional_boundary_theorem_only':True,
        'independent_model_reconstruction':True,'independent_fraction_intervals':True,
        'exact_binary_cover':True,'safe_denominators':True,
        'conditional_intersections_audited':True,'tree_nodes':len(nodes),
        'closed_leaves_verified':closed,'open_leaves':opened,'maximum_depth':maxdepth,
        'recorded_reason_matches':dict(recorded_reasons)==dict(independent_reasons),
        'independent_reasons':dict(independent_reasons),'seconds':time.monotonic()-start,
        'not_a_full_preclip_projection_theorem':True}
    target=Path(__file__).with_name('preclip-below-independent-replay.json')
    target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
