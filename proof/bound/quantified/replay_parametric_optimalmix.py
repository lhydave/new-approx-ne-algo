#!/usr/bin/env python3
"""Independent exact checker of the generic OM compilation and old fixture.

Does not import the compiler or any legacy exclusion evaluator.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,json
import sympy as S

HERE=Path(__file__).resolve().parent


def run():
    report=json.loads((HERE/'parametric-optimalmix-compilation.json').read_text())
    rho,sigma,eta,xi,g,r,lambda_payoff,mu_payoff,p,q,theta=S.symbols('rho sigma eta xi g r lambda_payoff mu_payoff p q theta')
    symbols={str(z):z for z in [rho,sigma,eta,xi,g,r,lambda_payoff,mu_payoff,p,q,theta]}
    def parse(x):return S.sympify(x,locals=symbols)
    checked=[]
    def equal(name,l,r):
        assert S.expand(S.together(l-r).as_numer_denom()[0])==0,name
        checked.append(name)
    aa,bb,t=S.symbols('aa bb t',positive=True)
    cut=bb/(aa+bb)
    equal('generic_ratio_intersection',(1-cut)/aa,cut/bb)
    equal('generic_ratio_order_factor',aa*bb*((1-t)/aa-t/bb),(aa+bb)*(cut-t))
    # This sign factor proves selected branches on both closed subintervals;
    # l is unrestricted, so the signed intercept is multiplied only afterward.
    generic_l,A,B,k=S.symbols('signed_l A B k')
    equal('generic_nonnegative_remainder',(1-t)*A+t*B-k*(aa*A+bb*B),
          ((1-t)-k*aa)*A+(t-k*bb)*B)
    p_dagger=parse(report['automatically_solved_cuts']['p_dagger'])
    q_dagger=parse(report['automatically_solved_cuts']['q_dagger'])
    equal('p_cut_satisfies_generic_ratio_equation',(1-p_dagger)/(eta*sigma/rho),p_dagger/(1-eta))
    equal('q_cut_satisfies_generic_ratio_equation',(1-q_dagger)/(xi*rho/sigma),q_dagger/(1-xi))
    # Denominators are positive under rho,sigma>0, eta,xi in (0,1).
    equal('p_cut_normalized_denominator',p_dagger,rho*(1-eta)/(rho*(1-eta)+sigma*eta))
    equal('q_cut_normalized_denominator',q_dagger,sigma*(1-xi)/(sigma*(1-xi)+rho*xi))
    lc=r-eta*(1-g/rho);lr=r-xi*(1-g/sigma)
    ps=[S.Integer(0),p_dagger,S.Integer(1)];qs=[S.Integer(0),q_dagger,S.Integer(1)]
    assert len(report['edges'])==12
    actualids=set()
    for edge in report['edges']:
        orient=edge['orientation'];level=edge['level'];index=edge['interval']
        eid=orient[0]+'_'+str(level)+'_'+str(index)
        assert edge['id']==eid and eid not in actualids
        actualids.add(eid)
        if orient=='vertical':
            pp=ps[level];qq=qs[index]+theta*(qs[index+1]-qs[index])
            i=min(level,1);j=index
        else:
            pp=ps[index]+theta*(ps[index+1]-ps[index]);qq=qs[level]
            i=index;j=min(level,1)
        wc=p/(1-eta) if i==0 else (1-p)*rho/(eta*sigma)
        wr=q/(1-xi) if j==0 else (1-q)*sigma/(xi*rho)
        R=(1-p)*(1-q)*g+q*(1-p*lambda_payoff)-lr*wr
        C=(1-p)*(1-q)*g+p*(1-q*mu_payoff)-lc*wc
        rr=S.cancel(R.subs({p:pp,q:qq}));cc=S.cancel(C.subs({p:pp,q:qq}))
        for expression in [rr,cc]:
            assert S.Poly(S.together(expression).as_numer_denom()[0],theta).degree()<=1
        vals={key:parse(val) for key,val in edge['caps'].items()}
        for key,value in [('a_R',rr.subs(theta,0)),('b_R',rr.subs(theta,1)),
                          ('a_C',cc.subs(theta,0)),('b_C',cc.subs(theta,1))]:
            equal(eid+'_'+key,vals[key],value)
        equal(eid+'_whole_R_segment',rr,(1-theta)*vals['a_R']+theta*vals['b_R'])
        equal(eid+'_whole_C_segment',cc,(1-theta)*vals['a_C']+theta*vals['b_C'])
        equal(eid+'_D',parse(edge['D']),(vals['b_R']-vals['a_R'])-(vals['b_C']-vals['a_C']))
        equal(eid+'_N',parse(edge['N']),vals['a_C']-vals['a_R'])
        equal(eid+'_W',parse(edge['W']),vals['b_R']*vals['a_C']-vals['a_R']*vals['b_C'])
    assert actualids=={o+'_'+str(i)+'_'+str(j)for o in ['v','h']for i in range(3)for j in range(2)}

    # Compare automatically generated internal edges to adopted old scalar caps.
    D_c=rho*(1-eta)+sigma*eta;D_r=sigma*(1-xi)+rho*xi
    J_C=rho*(r-eta*(1-r/rho))/D_c;J_R=sigma*(r-xi*(1-r/sigma))/D_r
    knot=((1-p_dagger)*(1-q_dagger)*r+q_dagger*(1-p_dagger*r/rho)-J_R,
          (1-p_dagger)*(1-q_dagger)*r+p_dagger*(1-q_dagger*r/sigma)-J_C)
    endpoints={
        'v_1_0':(((1-p_dagger)*r,(1-p_dagger)*r+p_dagger-J_C),knot),
        'v_1_1':(knot,(1-p_dagger*r/rho,p_dagger*(1-r/sigma)-J_C)),
        'h_1_0':(((1-q_dagger)*r+q_dagger-J_R,(1-q_dagger)*r),knot),
        'h_1_1':(knot,(q_dagger*(1-r/rho)-J_R,1-q_dagger*r/sigma))}
    for edge in report['edges']:
        if edge['id']not in endpoints:continue
        (R0,C0),(R1,C1)=endpoints[edge['id']]
        for key,value in [('a_R',R0),('b_R',R1),('a_C',C0),('b_C',C1)]:
            equal('adopted_'+edge['id']+'_'+key,
                  parse(edge['caps'][key]).subs({g:r,lambda_payoff:r/rho,mu_payoff:r/sigma}),value)

    # Reuse one EXISTING exact fixture, no new search or optimization.
    fixture_path=HERE.parent/'critical/alignment-exact-fixture.json'
    original={k:F(v)for k,v in json.loads(fixture_path.read_text())['coordinates'].items()}
    def segment_min(v0,v1):
        R0,C0=v0;R1,C1=v1
        values=[max(v0),max(v1)];den=(R1-R0)-(C1-C0)
        if den:
            t=(C0-R0)/den
            if 0<=t<=1:values.append((1-t)*R0+t*R1)
        return min(values)
    def verify_fixture(t):
        r=t['r'];rho=t['rho'];sigma=1-rho;eta=t['eta'];xi=t['xi'];a=t['alpha'];b=t['gamma']
        L=r/rho;M=r/sigma
        vr=min(F(1),(1-eta-r)/eta);vc=min(F(1),(1-xi-r)/xi)
        s_1=min(F(1),1-L-(r-(1-eta)*(1-M))/eta)
        s_2=min(F(1),1-M-(r-(1-xi)*(1-L))/xi)
        kappa_1=max(F(0),L-sigma/rho*vc);kappa_2=max(F(0),M-rho/sigma*vr)
        h_1=1-kappa_1;h_2=1-kappa_2
        capS=min(F(1),h_1-(r-(1-eta)*s_2)/eta)
        capT=min(F(1),h_2-(r-(1-xi)*s_1)/xi)
        domains=[rho-r,sigma-r,eta-r,1-r-eta,xi-r,1-r-xi,a-r,1-r-a,b-r,1-r-b]
        nonnegative=[vr,vc,s_1,s_2,capS,capT,
           eta*sigma/rho-(r-eta*(1-L)),xi*rho/sigma-(r-xi*(1-M))]
        firstrow=[(F(0),F(1)),(1-L,1-M),(h_1,s_2)]
        firstcol=[(F(1),F(0)),(1-L,1-M),(s_1,h_2)]
        secondrow=[(vr,F(1)),(s_1,h_2),(capS,capT),(F(1),F(0))]
        secondcol=[(F(1),vc),(h_1,s_2),(capS,capT),(F(0),F(1))]
        source=[eta*R+(1-eta)*C-r for R,C in firstrow]\
             +[(1-xi)*R+xi*C-r for R,C in firstcol]\
             +[a*R+(1-a)*C-r for R,C in secondrow]\
             +[(1-b)*R+b*C-r for R,C in secondcol]
        v=(1-xi)/xi;u=(1-eta)/eta
        NR=a*((1-b)*capS+b*(1-r/xi+v*s_1)-r)+(1-a)*((1-b)*capT+b*v*h_2)-r*(1-b+b*v)
        NC=(1-b)*((1-a)*capS+a*u*h_1)+b*((1-a)*capT+a*(1-r/eta+u*s_2)-r)-r*(1-a+a*u)
        corners=[(r,r),(F(0),F(1)),(F(1),F(0)),(1-L,1-M)]
        outer=[segment_min(corners[i],corners[j])-r for i,j in [(0,1),(0,2),(1,3),(2,3)]]
        assert all(z>=0 for z in domains+nonnegative+source+[NR,NC]+outer)
        return {'r':str(r),'rho':str(rho),'first4_and_shared_envelopes_pass':True,
                'NC_original_outer_chord_residuals':[str(z)for z in outer],
                'row_envelope':str(NR),'column_envelope':str(NC)}
    existing=verify_fixture(original)
    mirrored={**original,'rho':1-original['rho'],'eta':original['xi'],'xi':original['eta'],
              'alpha':original['gamma'],'gamma':original['alpha']}
    canonical=verify_fixture(mirrored)
    assert mirrored['rho']<=F(1,2)
    return {'success':True,'symbolic_identity_count':len(checked),'identity_ids':checked,
        'cells':4,'edges':12,'all_generated_segments_affine':True,
        'automatically_derived_old_internal_edges':4,
        'signed_intercept_requires_no_sign_assumption':True,
        'fixture_source':str(fixture_path),'search_used':False,
        'original_outer_edge_only_relaxation':{'existing':existing,'canonical_player_swap':canonical,
            'conclusion':'At r=15477/50000>r_star the adopted cap relaxation without internal edges remains feasible; original outer-edge-only encoding cannot close this backend',
            'not_claimed':'No claim that this scalar fixture is a game or that every stronger full shared-payoff NC encoding must fail'},
        'scope':'Generic linearization/cut/cell/edge identities and an existing exact fixture; not game realizability or full joint mass QE'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path)
    args=parser.parse_args();result=run()
    if args.out:args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'success':True,'symbolic_identity_count':result['symbolic_identity_count'],
                      'edges':12,'existing_fixture_passed':True,'out':str(args.out)if args.out else None}))
