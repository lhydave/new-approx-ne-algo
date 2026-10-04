#!/usr/bin/env python3
"""Small exact replay for the quantified/necessary interface.

No producer imports, CAD, optimization, random search or interval trees.
The checks prove displayed rational identities and exact finite fixtures;
they do not formalize LP duality, universal quantifier scope or convexity.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import sympy as sp

HERE = Path(__file__).resolve().parent


def run():
    checks = []

    def identity(name, left, right, assumptions=()):
        numerator = sp.together(left-right).as_numer_denom()[0]
        assert sp.expand(numerator) == 0, (name, sp.factor(numerator))
        checks.append({"id": name, "kind": "exact rational-function identity",
                       "equal_after_positive_denominator_clearing": True,
                       "domain_assumptions": list(assumptions)})

    rho,sigma,eta,xi,alpha,gamma,r,u,v = sp.symbols(
        'rho sigma eta xi alpha gamma r u v')
    R_x_s_1_x_s_2,R_x_x_s_2,R_y_s_1_y,R_x_s_1_y,C_x_y_s_2,C_x_x_s_2,C_x_s_1_y,C_x_s_1_x_s_2 = sp.symbols(
        'R_x_s_1_x_s_2 R_x_x_s_2 R_y_s_1_y R_x_s_1_y C_x_y_s_2 C_x_x_s_2 C_x_s_1_y C_x_s_1_x_s_2')
    A_s_comparison = -rho*R_x_x_s_2 + sigma*(C_x_y_s_2-C_x_x_s_2)
    B_s_comparison = -sigma*C_x_s_1_y + rho*(R_y_s_1_y-R_x_s_1_y)
    joint = rho*(R_y_s_1_y-R_x_s_1_y-R_x_x_s_2+R_x_s_1_x_s_2)+sigma*(C_x_y_s_2-C_x_x_s_2-C_x_s_1_y+C_x_s_1_x_s_2)
    identity('T04_T05_raw_joint_plane', rho*R_x_s_1_x_s_2+sigma*C_x_s_1_x_s_2+A_s_comparison+B_s_comparison, joint)

    h,f_1_value,f_2_value = sp.symbols('h f_1_value f_2_value')
    identity('R01_R02_nonnegative_slack_sum',
             h-(eta*f_1_value+(1-eta)*f_2_value), eta*(h-f_1_value)+(1-eta)*(h-f_2_value))
    identity('C01_C02_nonnegative_slack_sum',
             h-((1-xi)*f_1_value+xi*f_2_value), (1-xi)*(h-f_1_value)+xi*(h-f_2_value))

    A,D = sp.symbols('A D')
    D_c = rho*(1-eta)+sigma*eta
    p_dagger = rho*(1-eta)/D_c
    identity('proof_a_balances_same_column_deficit_plane',
             (1-p_dagger)*A+p_dagger*D, rho/D_c*(eta*(sigma/rho)*A+(1-eta)*D),
             ('rho>0,sigma>0,0<=eta<=1',))
    D_r = sigma*(1-xi)+rho*xi
    q_dagger = sigma*(1-xi)/D_r
    identity('proof_b_balances_same_row_deficit_plane',
             (1-q_dagger)*A+q_dagger*D, sigma/D_r*(xi*(rho/sigma)*A+(1-xi)*D),
             ('rho>0,sigma>0,0<=xi<=1',))

    S,T,s_1,s_2,h_1,h_2 = sp.symbols('S T s_1 s_2 h_1 h_2')
    H=alpha*S+(1-alpha)*T
    U=alpha*s_1+(1-alpha)*h_2
    K=(1-gamma)*S+gamma*T
    V=(1-gamma)*h_1+gamma*s_2
    denR=1-gamma+gamma*v
    denC=1-alpha+alpha*u
    E_1=alpha*((1-gamma)*S+gamma*(1-r*(1+v)+v*s_1)-r)\
        +(1-alpha)*((1-gamma)*T+gamma*v*h_2)-r*denR
    E_2=(1-gamma)*((1-alpha)*S+alpha*u*h_1)\
        +gamma*((1-alpha)*T+alpha*(1-r*(1+u)+u*s_2)-r)-r*denC
    surplusR=(1-gamma)*(H-r)+gamma*v*(U-r)\
        +alpha*gamma*(1-r-r*v)-alpha*r
    surplusC=(1-alpha)*(K-r)+alpha*u*(V-r)\
        +alpha*gamma*(1-r-r*u)-gamma*r
    identity('same_alpha_gamma_row_surplus',E_1,surplusR)
    identity('same_alpha_gamma_column_surplus',E_2,surplusC)
    assert sp.Poly(E_1,alpha,gamma).degree(alpha)<=1
    assert sp.Poly(E_1,alpha,gamma).degree(gamma)<=1
    assert sp.Poly(E_2,alpha,gamma).degree(alpha)<=1
    assert sp.Poly(E_2,alpha,gamma).degree(gamma)<=1
    checks.append({"id":"cleared_envelopes_are_joint_multiaffine",
                   "kind":"symbolic degree check","passed":True})

    al,ah,gl,gh,s,t,c00,c10,c01,c11 = sp.symbols(
        'al ah gl gh s t c00 c10 c01 c11')
    def bilinear(x,y):
        return c00+c10*x+c01*y+c11*x*y
    barycentric = (1-s)*(1-t)*bilinear(al,gl)+s*(1-t)*bilinear(ah,gl)\
                  +(1-s)*t*bilinear(al,gh)+s*t*bilinear(ah,gh)
    identity('complete_bilinear_corner_interpolation',
             bilinear(al+s*(ah-al),gl+t*(gh-gl)),barycentric,
             ('0<=s,t<=1 gives nonnegative corner weights',))

    a_R,b_R,a_C,b_C = sp.symbols('a_R b_R a_C b_C')
    d=(b_R-a_R)-(b_C-a_C)
    n=a_C-a_R
    w=b_R*a_C-a_R*b_C
    theta=n/d
    identity('segment_crossing_row_is_W_over_D', a_R+(b_R-a_R)*theta, w/d,
             ('D!=0','crossing guarded to theta in [0,1]'))
    identity('segment_crossing_column_is_W_over_D', a_C+(b_C-a_C)*theta, w/d,
             ('D!=0','crossing guarded to theta in [0,1]'))

    def exact_minimum(ar,br,ac,bc):
        # Independently evaluate endpoints plus actual line intersection.
        values=[max(ar,ac),max(br,bc)]
        slope_difference=(br-ar)-(bc-ac)
        if slope_difference:
            crossing=(ac-ar)/slope_difference
            if 0<=crossing<=1:
                values.append(max(ar+(br-ar)*crossing,ac+(bc-ac)*crossing))
        return min(values)

    def eliminated_formula(ar,br,ac,bc,threshold):
        d=(br-ar)-(bc-ac);n=ac-ar;w=br*ac-ar*bc
        if max(ar,ac)<threshold or max(br,bc)<threshold:
            return False
        if d>0 and 0<=n<=d:
            return w-threshold*d>=0
        if d<0 and d<=n<=0:
            return w-threshold*d<=0
        return True

    fixtures=0
    branches={'parallel':0,'positive_crossing':0,'negative_crossing':0,'outside_crossing':0}
    endpoints=[F(-1),F(0),F(1,2),F(1)]
    thresholds=[F(-1,3),F(0),F(1,3),F(1)]
    for ar in endpoints:
        for br in endpoints:
            for ac in endpoints:
                for bc in endpoints:
                    d=(br-ar)-(bc-ac);n=ac-ar
                    key='parallel' if d==0 else\
                        'positive_crossing' if d>0 and 0<=n<=d else\
                        'negative_crossing' if d<0 and d<=n<=0 else 'outside_crossing'
                    for threshold in thresholds:
                        assert eliminated_formula(ar,br,ac,bc,threshold)==\
                            (exact_minimum(ar,br,ac,bc)>=threshold)
                        branches[key]+=1;fixtures+=1
    checks.append({"id":"complete_segment_quantifier_cases",
                   "kind":"exact finite rational fixtures",
                   "fixtures":fixtures,"branches":branches,
                   "not_a_substitute_for_universal_piecewise_affine_lemma":True})

    # Explicitly demonstrates why individual corner maxima are not joint QE.
    cornerR=max(g-F(3,5) for g in [F(0),F(1)])
    cornerC=max(F(2,5)-g for g in [F(0),F(1)])
    assert cornerR>=0 and cornerC>=0 and F(3,5)>F(2,5)
    checks.append({"id":"separate_maxima_not_joint_mass_feasibility",
                   "kind":"exact rational counterexample",
                   "toy_row_surplus":"gamma-3/5","toy_column_surplus":"2/5-gamma",
                   "separate_maxima":[str(cornerR),str(cornerC)],
                   "joint_lower_gamma":"3/5","joint_upper_gamma":"2/5",
                   "joint_interval_empty":True})

    # Checks the interface document's Stage-I boundary; not mathematical proof.
    registry=json.loads((HERE/'quantified-registry.json').read_text())
    assert registry['stage_i']['target_parameters']==[]
    assert registry['stage_i']['explicit_special_mixing_assignments']==[]
    assert registry['stage_i']['early_branches']==[]
    assert registry['stage_i']['primitive_blocks']==[
        'TSenlarged','FixedRowMix','FixedColumnMix','OptimalMixing','SelectBest']
    assert registry['stage_i']['ssa_order'][-1]=='o=SelectBest(first4,O)'
    assert registry['stage_i']['candidate_profiles']==[
        ['y_s_1','x_1_2'],['x_1_1','y_s_2'],['y_1_1','x_2_2'],['x_2_1','y_1_2'],['x_o','y_o']]
    clauses={clause['id']:clause for clause in registry['stage_i']['clauses']}
    assert clauses['O02']['property']=='F(x_o,y_o)<=min(F(x_s_1,x_s_2),F(x_s_1,y_s_2),F(y_s_1,x_s_2),F(y_s_1,y_s_2))'
    assert clauses['B01']['property']=='output belongs to five actual profiles [first4,O] and has F no greater than every list member'
    assert registry['stage_ii']['exact_joint_mass_QE']['implemented_here'] is False
    contracts=json.loads((HERE/'stage-i-contracts.json').read_text())
    assert contracts['schema']==1
    assert contracts['scope']=='Stage-I building-block contracts and SSA composition'
    first=contracts['contract_text']
    for forbidden in ['ρ','c+δ','p_dagger=rho','q_dagger=sigma','T_a=','T_b=','D_c=']:
        assert forbidden not in first, forbidden
    checks.append({"id":"stage_I_contains_only_final_block_contracts",
                   "kind":"registry/text boundary check","passed":True})
    return {"success":True,"checks":checks,"check_count":len(checks),
            "scope":"Only rational identities, complete finite segment fixtures and registry boundaries",
            "does_not_certify":["LP strong-duality theorem","general convexity theorem",
                 "all original-game quantifier discharge","game realizability",
                 "complete joint alpha/gamma QE","natural bound global certificates"],
            "uses_floating_point":False,"runs_CAD":False,"replays_large_trees":False}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    result=run()
    output=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.out:
        args.out.write_text(output)
    print(json.dumps({"success":result['success'],"check_count":result['check_count'],
                      "segment_fixtures":1024,"out":str(args.out) if args.out else None}))
