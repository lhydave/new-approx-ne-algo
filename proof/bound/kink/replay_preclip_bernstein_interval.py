"""Independent replay of global/preclip-bernstein-certificate.json.

Uses separate Fraction intervals, independent symbolic arc clearing, direct
power-to-Bernstein conversion, and closed-form half-interval transforms.
It never imports or calls the producer, and never adds tree nodes.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json,sys,time,hashlib
import sympy as s
from replay_preclip_interval import B

r,rho,u,v=s.symbols('r rho u v');sigma=1-rho;c=1-2*r
X=rho*(1-r)-r;Y=sigma*(1-r)-r
sr_num=sigma*X+rho*Y*u;sc_num=rho*Y+sigma*X*v
hc_num=sigma-r*(1+rho)+rho*(1-r)*u
hr_num=rho-r*(1+sigma)+sigma*(1-r)*v
M=sigma*hr_num-2*r*rho*sigma+u*(sc_num-r*rho*sigma)
N=2*r*rho*sigma-rho*hc_num-v*(sr_num-r*rho*sigma)
N_R=sigma*Y+rho*(sigma-2*r)*u+sigma*sigma*c*v+sigma*(X-r*rho)*u*v
N_C=rho*X+sigma*(rho-2*r)*v+rho*rho*c*u+rho*(Y-r*sigma)*u*v
L_scaled=(X+rho*c*u)*N_R+r*rho*u*N_C
L_scaled_0=L_scaled.subs(v,0);v_den=-s.diff(L_scaled,v)
def clear(expr):
    assert s.Poly(s.expand(expr),v).degree()<=1
    return s.expand(v_den*expr.subs(v,0)+L_scaled_0*s.diff(expr,v))
MA=clear(M);NA=clear(N)
polys={'v_den':v_den,'L_scaled_0':L_scaled_0,'MA':MA,'NA':NA,'UA':clear(rho*M*(hc_num-r*sigma)-N*(r*rho*sigma-sr_num)),
       'MN':MA-NA,'TA':r*rho*sigma*v_den-NA,'KUA':clear(s.diff(L_scaled,u)),
       'SCA':clear(sc_num),'SR':sr_num,'mass':rho-r,'umass':u*(1-r)-r,
       'zlower':(1-r)*NA-r*MA,'zupper':(1-r)*MA-r*NA,
       'clipgraph':sigma*(1-r)*L_scaled_0-r*(1+sigma)*v_den,'vsource':L_scaled_0-s.Rational(4,5)*v_den}
AA=hc_num-r*sigma;BB=2*r*rho*sigma-sr_num
EA=rho*L_scaled_0*AA*MA**2-NA*MA*(L_scaled_0*BB-c*rho*sigma*v_den)-r*rho*sigma*NA**2*v_den
polys['EA']=EA
# A generic identity plus independently checked linear arc clearing proves
# EA=EH*v_den^3 without relying on a large expanded graph numerator.
ap,aq,ar,ac,aa,bb,kk,ll,mm,nn=s.symbols('rho sigma R C A B L_scaled v_den M N')
generic_EH=ap*(kk/ll)*aa*(mm/ll)**2-(nn/ll)*(mm/ll)*(kk*bb/ll-ac*ap*aq)-ar*ap*aq*(nn/ll)**2
generic_EA=ap*kk*aa*mm**2-nn*mm*(kk*bb-ac*ap*aq*ll)-ar*ap*aq*nn**2*ll
assert s.cancel(generic_EH*ll**3-generic_EA)==0
assert s.expand(MA-v_den*M.subs(v,0)-L_scaled_0*s.diff(M,v))==0
assert s.expand(NA-v_den*N.subs(v,0)-L_scaled_0*s.diff(N,v))==0
# Independent exact algebra audit of the cleared T cut and tangent numerator.
assert s.simplify((r*rho*sigma-N)-rho*sigma*(hc_num/sigma-r+v*(sr_num/(rho*sigma)-r)))==0


def initial_grids(root):
    a=s.symbols('a:3');grids={};degrees={}
    mapping={var:s.Rational(lo)+(s.Rational(hi)-s.Rational(lo))*z for var,z,(lo,hi)in zip((r,rho,u),a,root)}
    for name,expr in polys.items():
        pol=s.Poly(s.expand(expr.subs(mapping)),*a);ds=tuple(pol.degree_list());power={idx:F(co)for idx,co in pol.terms()}
        # The tensor transform is applied in order 2,0,1; the independently
        # derived direct formula is checked at four distinct grid positions.
        bern=power
        for axis in (2,0,1):
            transformed={}
            for idx in product(*(range(d+1)for d in ds)):
                value=F(0)
                for j in range(idx[axis]+1):
                    key=list(idx);key[axis]=j
                    value+=bern.get(tuple(key),F(0))*F(comb(idx[axis],j),comb(ds[axis],j))
                transformed[idx]=value
            bern=transformed
        audit_indices={(0,0,0),ds,tuple(d//2 for d in ds),(ds[0],0,ds[2])}
        for idx in audit_indices:
            direct=F(0)
            for powers,co in power.items():
                if all(k<=i for k,i in zip(powers,idx)):
                    weight=F(1)
                    for k,i,d in zip(powers,idx,ds):weight*=F(comb(i,k),comb(d,k))
                    direct+=co*weight
            assert direct==bern[idx]
        grids[name]=bern;degrees[name]=ds
    return grids,degrees

HALF={}
def half_grid(grid,ds,axis):
    n=ds[axis]
    if n not in HALF:
        HALF[n]=([[F(comb(k,j),2**k)for j in range(k+1)]for k in range(n+1)],
                 [[F(comb(n-k,j-k),2**(n-k))for j in range(k,n+1)]for k in range(n+1)])
    wl,wr=HALF[n];others=[z for z in range(3)if z!=axis];left={};right={}
    for fixed in product(*(range(ds[z]+1)for z in others)):
        idx=[0]*3
        for z,i in zip(others,fixed):idx[z]=i
        row=[]
        for k in range(n+1):idx[axis]=k;row.append(grid[tuple(idx)])
        for k in range(n+1):
            idx[axis]=k
            left[tuple(idx)]=sum((row[j]*wl[k][j]for j in range(k+1)),F(0))
            right[tuple(idx)]=sum((row[j]*wr[k][j-k]for j in range(k,n+1)),F(0))
    return left,right

def intersection(a,b):
    lo=max(a.a,b.a);hi=min(a.b,b.b)
    return None if lo>hi else B(lo,hi)
def contract(a,lo=None,hi=None):
    low=a.a if lo is None else max(a.a,lo);high=a.b if hi is None else min(a.b,hi)
    return None if low>high else B(low,high)

def replay_leaf(node,box,grids,EA_bound):
    if node['kind']=='unresolved':return 'unresolved'
    bound={name:B(min(grid.values()),max(grid.values()))for name,grid in grids.items()}
    name=node['reason'];rr,pp,uu=[B(*z)for z in box];qq=1-pp
    v0=rr*(1+qq)/(qq*(1-rr))
    if name in ('bern_mass','bern_umass','bern_SR'):
        assert bound[name[5:]].b<0;return 'excluded'
    if name=='bern_L_scaled_domain':
        test=bound['L_scaled_0']-bound['v_den']*B(F(4,5),v0.b)
        assert test.a>0 or test.b<0;return 'excluded'
    # All homogeneous divisions, source/E1 cuts and contractions are guarded.
    assert bound['v_den'].a>0,('unexpected unguarded producer leaf',name)
    if name=='bern_EA':assert EA_bound().b<0;return 'excluded'
    if name=='bern_clipgraph':assert bound['clipgraph'].a>0;return 'excluded'
    if name=='bern_vsource':assert bound['vsource'].b<=0;return 'excluded'
    strict_names={'MA','NA','zlower','MN'}
    if name.startswith('bern_')and name[5:]in ('MA','NA','UA','TA','SCA','zupper','zlower','MN'):
        f=bound[name[5:]];assert f.b<0 or(name[5:]in strict_names and f.b<=0);return 'excluded'
    if name=='bern_upper_root':assert bound['KUA'].a>0;return 'excluded'
    r,rho,u=rr,pp,uu;sigma=1-rho;c=1-2*r;k=r/(1-r)
    vv=intersection(bound['L_scaled_0']/bound['v_den'],B(max(k.a,F(4,5)),v0.b))
    if name=='bern_v_domain':assert vv is None;return 'excluded'
    assert vv is not None;v=vv
    X=rho*(1-r)-r;Y=sigma*(1-r)-r
    A1=intersection(sigma*X+rho*Y*u,bound['SR'])
    H=sigma-r*(1+rho)+rho*(1-r)*u;Rh=rho-r*(1+sigma)+sigma*(1-r)*v
    B1=intersection(rho*Y+sigma*X*v,bound['SCA']/bound['v_den'])
    if A1 is None or B1 is None:
        assert name=='bern_intersection';return 'excluded'
    M=intersection(sigma*Rh-2*r*rho*sigma+u*(B1-r*rho*sigma),bound['MA']/bound['v_den'])
    N=intersection(2*r*rho*sigma-rho*H-v*(A1-r*rho*sigma),bound['NA']/bound['v_den'])
    if M is None or N is None:
        assert name=='bern_intersection';return 'excluded'
    # Derived M>N uses E>=0, U>=0 and the already-proved v_den>0 guard.
    # N<=rPQ is exactly T>=0, rather than the weaker N<=2rPQ.
    M=contract(M,lo=max(F(0),N.a))
    N=contract(N,lo=F(0),hi=min((r*rho*sigma).b,M.b))if M is not None else None
    if M is None or N is None:
        assert name=='bern_MN_T_contract';return 'excluded'
    SV=sigma*(1-r)+X*u;TU=rho*(1-r)+Y*v
    CS=r*rho*sigma-B1;DS=r*rho*sigma-A1;AA=H-r*sigma;BB=2*r*rho*sigma-A1;LL=v*BB-c*rho*sigma
    KU=bound['KUA']/bound['v_den'];KV=-bound['v_den']
    raw_EH=rho*v*AA*M**2-N*M*LL-r*rho*sigma*N**2
    if name=='E_nonnegative' and raw_EH.b<0:return 'excluded'
    strict=KU.b<0 and KV.b<0 and bound['umass'].a>0 and A1.a>0 and B1.a>0
    Z10=-M*v*BB-v*LL*sigma*SV-2*r*rho*sigma*v*DS
    Z11=-rho*M*v*v*Y-v*LL*CS-2*r*rho*rho*sigma*v*TU
    Z00=rho*M*v*AA+2*rho*v*v*AA*sigma*SV-v*LL*DS
    Z01=-rho*rho*M*v*v*(1-r)+2*rho*v*v*AA*CS-v*LL*rho*TU
    J1=(Z01*M+Z11*N)*M
    def target(eh):return -KU*((Z00*M+Z10*N)*M-2*rho*sigma*eh)-KV*J1
    if node['kind']=='proved':
        assert strict
        if name=='strict_slopes_caps_E_positive' and raw_EH.a>0:return 'proved'
        if name=='strict_slopes_caps_zero_derivative' and M.a>0 and target(raw_EH).a>0:return 'proved'
    # EA is propagated only when the terminal claim actually needs its bound.
    # This is a proof optimization, not prefix or report reuse.
    EH=intersection(raw_EH,EA_bound()/(bound['v_den']**3))
    if EH is None:
        assert node['kind']=='excluded';return 'excluded'
    if name=='E_nonnegative':assert EH.b<0;return 'excluded'
    assert node['kind']=='proved' and strict
    if name=='strict_slopes_caps_E_positive':assert EH.a>0;return 'proved'
    assert name=='strict_slopes_caps_zero_derivative'
    assert M.a>0 and target(EH).a>0
    return 'proved'

def main(path,deep=False):
    started=time.monotonic();raw=Path(path).read_bytes();doc=json.loads(raw);root=tuple(tuple(F(z)for z in pair)for pair in doc['root'])
    prefix_path=Path(__file__).with_name('preclip-bernstein-initial-4991.json')
    prefix_report_path=Path(__file__).with_name('preclip-bernstein-initial-4991-independent-replay.json')
    if deep:
        prefix_raw=None;prefix=None
    else:
        prefix_raw=prefix_path.read_bytes();prefix=json.loads(prefix_raw);report=json.loads(prefix_report_path.read_text())
        assert report['accepted_leaf_audit'] and report['coverage_partition_verified']
        assert report['certificate_sha256']==hashlib.sha256(prefix_raw).hexdigest()
        assert prefix['root']==doc['root']
    grids,degrees=initial_grids(root);assert {n:list(d)for n,d in degrees.items()}==doc['degrees']
    ea_root=grids.pop('EA');ea_cache={():ea_root};ea_audit={'requests':0,'bisections':0}
    def ea_grid(path):
        ea_audit['requests']+=1
        if path in ea_cache:return ea_cache[path]
        parent=path[:-1];axis,side=path[-1];vals=ea_grid(parent)
        left,right=half_grid(vals,degrees['EA'],axis);ea_audit['bisections']+=1
        ea_cache[parent+((axis,0),)]=left;ea_cache[parent+((axis,1),)]=right
        return ea_cache[path]
    stats={'excluded':0,'proved':0,'unresolved':0,'split':0};reasons={};last_report=time.monotonic();reused=0;fresh=0
    cache={}
    def old_metrics(n):
        if 'split'in n:
            ls,lr=old_metrics(n['left']);rs,rr=old_metrics(n['right']);ss={k:ls[k]+rs[k]for k in stats};ss['split']+=1;co=lr.copy()
            for k,z in rr.items():co[k]=co.get(k,0)+z
        else:
            ss={k:0 for k in stats};ss[n['kind']]=1;co={n.get('reason','unresolved'):1}
        cache[id(n)]=(ss,co);return ss,co
    if prefix is not None:old_metrics(prefix['tree'])
    def visit(node,box,gs,old=None,path=()):
        nonlocal last_report,reused,fresh
        if old is not None and cache[id(old)][0]['unresolved']==0:
            # Same subtree and same inherited bisection box: all terminal claims
            # were already independently replayed in the immutable prefix.
            assert node==old
            ss,co=cache[id(old)]
            for k,z in ss.items():stats[k]+=z
            for k,z in co.items():reasons[k]=reasons.get(k,0)+z
            reused+=ss['excluded']+ss['proved'];return
        if 'split'in node:
            axis=node['split'];assert axis in (0,1,2)
            if old is not None and 'split'in old:assert old['split']==axis
            a,b=box[axis];mid=(a+b)/2;left=list(box);right=list(box);left[axis]=(a,mid);right[axis]=(mid,b)
            gl={};gr={}
            for name,grid in gs.items():gl[name],gr[name]=half_grid(grid,degrees[name],axis)
            stats['split']+=1
            visit(node['left'],tuple(left),gl,old.get('left')if old is not None and 'split'in old else None,path+((axis,0),))
            visit(node['right'],tuple(right),gr,old.get('right')if old is not None and 'split'in old else None,path+((axis,1),))
        else:
            def get_EA():
                vv=ea_grid(path);return B(min(vv.values()),max(vv.values()))
            kind=replay_leaf(node,box,gs,get_EA);stats[kind]+=1;fresh+=1;key=node.get('reason','unresolved');reasons[key]=reasons.get(key,0)+1
            if time.monotonic()-last_report>30:
                print(json.dumps(dict(progress=stats,reused=reused,fresh=fresh,seconds=time.monotonic()-started)),flush=True);last_report=time.monotonic()
    visit(doc['tree'],root,grids,prefix['tree']if prefix is not None else None)
    assert stats['excluded']==doc['stats']['excluded'] and stats['proved']==doc['stats']['proved'] and stats['unresolved']==doc['unresolved']
    assert stats['split']+stats['excluded']+stats['proved']+stats['unresolved']==doc['stats']['nodes']
    result={'status':'complete_verified'if not stats['unresolved']else'partial_leaves_verified_no_global_theorem','certificate_sha256':hashlib.sha256(raw).hexdigest(),'accepted_leaf_audit':True,'coverage_partition_verified':True,'independent_symbolic_arc_clearing':True,'EA_equals_EH_times_Lcubed_symbolically_verified':True,'independent_Bernstein_conversion_and_half_transform':True,'deep_fresh_replay':deep,'EA_lazy_propagation':ea_audit,'cached_immutable_prefix_sha256':hashlib.sha256(prefix_raw).hexdigest()if prefix_raw is not None else None,'cached_prefix_accepted_leaves':reused,'new_terminal_claims_independently_checked':fresh,'upper_branch_assumption':'KU<=0 including KU=0 boundary','slope_and_cap_claim':'strict KU<0,KV<0,u>k,s_1>0,s_2>0 throughout the feasible E>=0 upper-knot chart','derivative_claim':'Jzero>0 only when E=0; offzero Jzero-2E is not the full derivative','guard_audit':'v_den>0 verified before all homogeneous cuts, z<1 contraction, quotient bounds or EA/v_den^3','T_cut_audit':'TA=rPQ*v_den-NA exactly encodes T>=0','stats':stats,'reasons':reasons,'elapsed_seconds':time.monotonic()-started}
    output=Path(path).with_name(Path(path).stem+('-deep-independent-replay.json'if deep else'-independent-replay.json'));output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main(sys.argv[1],deep='--deep'in sys.argv[2:])
