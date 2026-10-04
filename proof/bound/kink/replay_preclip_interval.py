"""Independent exact replay: SymPy-derived expressions, separate interval arithmetic.
It verifies every accepted leaf and reconstructs the complete binary partition.
An unresolved leaf prevents a global theorem, even if all accepted leaves pass.
"""
from fractions import Fraction as ExactFraction
from pathlib import Path
import sympy as s
import json,sys,time
from math import comb

class B:
    def __init__(self,a,b=None):self.a=ExactFraction(a);self.b=ExactFraction(a if b is None else b);assert self.a<=self.b
    @staticmethod
    def c(v):return v if isinstance(v,B)else B(v)
    def __neg__(self):return B(-self.b,-self.a)
    def __add__(self,v):v=B.c(v);return B(self.a+v.a,self.b+v.b)
    __radd__=__add__
    def __sub__(self,v):return self+-B.c(v)
    def __rsub__(self,v):return B.c(v)+-self
    def __mul__(self,v):
        v=B.c(v);rho=[self.a*v.a,self.a*v.b,self.b*v.a,self.b*v.b];return B(min(rho),max(rho))
    __rmul__=__mul__
    def __truediv__(self,v):
        v=B.c(v);assert v.b<0 or v.a>0
        return self*B(ExactFraction(1)/v.b,ExactFraction(1)/v.a)
    def __rtruediv__(self,v):return B.c(v)/self
    def __pow__(self,n):
        assert isinstance(n,int)
        if n==0:return B(1)
        if n<0:return B(1)/(self**(-n))
        if n%2==0:
            vals=[self.a**n,self.b**n];return B(0 if self.a<=0<=self.b else min(vals),max(vals))
        return B(self.a**n,self.b**n)

r,rho,u,v,t=s.symbols('r rho u v t');sigma=1-rho;c=1-2*r;k=r/(1-r)
X=rho*(1-r)-r;Y=sigma*(1-r)-r
A1=sigma*X+rho*Y*u;B1=rho*Y+sigma*X*v
N_R_0=sigma*Y+rho*(sigma-2*r)*u;N_R_1=sigma*sigma*c+sigma*(X-r*rho)*u
N_C_0=rho*X+rho*rho*c*u;N_C_1=sigma*(rho-2*r)+rho*(Y-r*sigma)*u
L_scaled_0=(X+rho*c*u)*N_R_0+r*rho*u*N_C_0
L_scaled_1=(X+rho*c*u)*N_R_1+r*rho*u*N_C_1
L_scaled=L_scaled_0+L_scaled_1*v;KU=s.diff(L_scaled,u);KV=s.diff(L_scaled,v)
H=sigma-r*(1+rho)+rho*(1-r)*u;R=rho-r*(1+sigma)+sigma*(1-r)*v
M=sigma*R-2*r*rho*sigma+u*(B1-r*rho*sigma);N=2*r*rho*sigma-rho*H-v*(A1-r*rho*sigma)
AA=H-r*sigma;DS=r*rho*sigma-A1;CS=r*rho*sigma-B1;BB=2*r*rho*sigma-A1;LL=v*BB-c*rho*sigma
SV=sigma*(1-r)+X*u;TU=rho*(1-r)+Y*v
EH=rho*v*AA*M**2-N*M*LL-r*rho*sigma*N**2
Z1=-M*v*BB-t*rho*M*v*v*Y-v*LL*(sigma*SV+t*CS)-2*r*rho*sigma*v*(DS+t*rho*TU)
Z0=rho*M*v*AA-t*rho*rho*M*v*v*(1-r)+2*rho*v*v*AA*(sigma*SV+t*CS)-v*LL*(DS+t*rho*TU)
ZH=(Z0*M+Z1*N)*M-2*rho*sigma*EH
TARGET=-KU*ZH.subs(t,0)-KV*s.diff(ZH,t)
F0=s.lambdify((r,rho,u),L_scaled_0,modules='math');F1=s.lambdify((r,rho,u),L_scaled_1,modules='math')
Fk=s.lambdify(r,k,modules='math');Fv0=s.lambdify((r,rho),r*(1+sigma)/(sigma*(1-r)),modules='math')
AC=[X*sigma*Y,rho*c*sigma*Y+X*rho*(sigma-2*r)+r*rho*rho*X,rho*rho*c*(sigma-2*r+r*rho)]
L_scaled_1_coefficients=[X*sigma*sigma*c,rho*c*sigma*sigma*c+X*sigma*(X-r*rho)+r*rho*sigma*(rho-2*r),rho*c*sigma*(X-r*rho)+r*rho*rho*(Y-r*sigma)]
assert all(s.expand(s.Poly(L_scaled_0,u).nth(i)-AC[i])==0 for i in range(3))
assert all(s.expand(s.Poly(L_scaled_1,u).nth(i)-L_scaled_1_coefficients[i])==0 for i in range(3))
ACF=[s.lambdify((r,rho),e,modules='math')for e in AC];BCF=[s.lambdify((r,rho),e,modules='math')for e in L_scaled_1_coefficients]

def bernstein(coeff,ub):
    n=len(coeff)-1;lo=ub.a;wid=ub.b-ub.a
    power=[sum((coeff[j]*comb(j,i)*lo**(j-i)*wid**i for j in range(i,n+1)),B(0))for i in range(n+1)]
    bez=[sum((power[i]*ExactFraction(comb(j,i),comb(n,i))for i in range(j+1)),B(0))for j in range(n+1)]
    return B(min(z.a for z in bez),max(z.b for z in bez))

expressions={'mass':rho-r,'u_mass':u-k,'s_1':A1,'v_mass':v-k,'v_clip':r*(1+sigma)/(sigma*(1-r))-v,'s_2':B1,'m':M,'n':N,'T':2*r*rho*sigma-N,'U':rho*M*AA-N*DS,'upper_root':-KU,'upper_root_closed':-KU,'lower_z':N-k*M,'upper_z_guarded_Kv':M-N,'E_zero':EH,'E_nonnegative':EH,'target':TARGET,'slope_u':-KU,'slope_v':-KV}
FF={name:s.lambdify((r,rho,u,v),expr,modules='math')for name,expr in expressions.items()}

def verify_leaf(node,box,modern=False):
    if node['kind']=='unresolved':return 'unresolved'
    rv,pv,uv=[B(*x)for x in box];kk=Fk(rv);v0=Fv0(rv,pv);kd=B(max(kk.a,ExactFraction(4,5))if modern else kk.a,v0.b)
    name=node['reason']
    if name in ('mass','u_mass','s_1'):
        bound=FF[name](rv,pv,uv,kd);assert bound.b<0;return 'excluded'
    if modern:
        aa=[f(rv,pv)for f in ACF];bb=[f(rv,pv)for f in BCF]
        k0=bernstein(aa,uv);k1=bernstein(bb,uv)
    else:k0=F0(rv,pv,uv);k1=F1(rv,pv,uv)
    if name in ('clip_K_boundary','source_K_boundary'):
        assert k1.b<0
        if name=='clip_K_boundary':
            cf=[(1-pv)*(1-rv)*aa[i]+rv*(2-pv)*bb[i]for i in range(3)]
            assert bernstein(cf,uv).a>0
        else:assert bernstein([aa[i]+ExactFraction(4,5)*bb[i]for i in range(3)],uv).b<=0
        return 'excluded' 
    if name=='K_domain':
        f=k0+k1*kd;assert f.b<0 or f.a>0;return 'excluded'
    assert k1.b<0 or k1.a>0
    vb=-k0/k1
    if name=='v_intersection':assert max(vb.a,kd.a)>min(vb.b,v0.b);return 'excluded'
    if name in ('v_mass','v_clip'):
        f=FF[name](rv,pv,uv,vb);assert f.b<0;return 'excluded'
    vb=B(max(vb.a,kd.a),min(vb.b,v0.b))
    if node['kind']=='excluded':
        f=FF[name](rv,pv,uv,vb)
        if name=='E_zero':assert f.b<0 or f.a>0
        elif name in ('m','n','upper_root','lower_z','upper_z_guarded_Kv'):assert f.b<=0
        else:assert f.b<0
        if name=='upper_z_guarded_Kv':assert k1.b<0
        return 'excluded'
    assert FF['slope_u'](rv,pv,uv,vb).a>0
    assert FF['slope_v'](rv,pv,uv,vb).a>0
    if modern:
        assert FF['u_mass'](rv,pv,uv,vb).a>0
        assert FF['s_1'](rv,pv,uv,vb).a>0
        assert FF['s_2'](rv,pv,uv,vb).a>0
    if name=='strict_slopes_caps_E_positive':assert FF['E_zero'](rv,pv,uv,vb).a>0
    else:
        assert FF['m'](rv,pv,uv,vb).a>0
        assert FF['target'](rv,pv,uv,vb).a>0
    return 'proved'

def main(path):
    start=time.monotonic();doc=json.loads(Path(path).read_text());root=tuple(tuple(ExactFraction(v)for v in z)for z in doc['root']);stats={'excluded':0,'proved':0,'unresolved':0,'split':0}
    def rec(n,box):
        if 'split'in n:
            d=n['split'];assert d in (0,1,2);a,b=box[d];mid=(a+b)/2;left=list(box);right=list(box);left[d]=(a,mid);right[d]=(mid,b);stats['split']+=1;rec(n['left'],tuple(left));rec(n['right'],tuple(right))
        else:stats[verify_leaf(n,box,'cuts'in doc)]+=1
    rec(doc['tree'],root)
    result=dict(status='complete_verified'if not stats['unresolved']else'partial_leaves_verified_no_global_proof',coverage_partition_verified=True,independent_symbolic_derivatives=True,stats=stats,elapsed_seconds=time.monotonic()-start)
    Path(path).with_name(Path(path).stem+'-replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main(sys.argv[1])
