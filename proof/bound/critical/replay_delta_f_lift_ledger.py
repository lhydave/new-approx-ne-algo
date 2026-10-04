"""Independent exact replay of the finite retained near-BR delta_f ledger.

No producer import, optimizations, root computations, or tree replays.
This proves arithmetic implications given the explicitly typed premises;
it does not formalize game semantics or contract implementation refinement.
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import re
import sympy as s

base=Path(__file__).parent
path=base/'delta_f-lift-ledger.json'
cert=json.loads(path.read_text())
registry=json.loads((base/cert['registry']).read_text())
valid_ids={x['id'] for x in registry['stage_i']['clauses']}
names=set()
for collection in [cert['nonnegative_atoms'],cert['premises']]:
 for item in collection.values():
  names.update(re.findall(r'[A-Za-z][A-Za-z0-9_]*',item['expression']))
for item in cert['nonnegative_combinations']:
 names.update(re.findall(r'[A-Za-z][A-Za-z0-9_]*',item['target']))
for expr in list(cert['zero_identities'].values())+list(cert['delta_f_free_ideal_aliases'].values())+list(cert['algebraic_definitions'].values()):
 names.update(re.findall(r'[A-Za-z][A-Za-z0-9_]*',expr))
names.update(cert['algebraic_definitions'])
symbols={name:s.Symbol(name,real=True) for name in names}
symbols.update(Min=s.Min,Max=s.Max)
definitions={}
def parse(expr):return s.sympify(expr,locals=symbols).xreplace(definitions)
def equal(a,b):assert s.factor(s.together(a-b))==0
expected_definitions={
 'g':'G_s+delta_f','R_x_s_1_x_s_2':'UR_x_s_2-g','R_y_s_1_x_s_2':'UR_x_s_2-tau_1',
 'C_x_s_1_x_s_2':'UC_x_s_1-g','C_x_s_1_y_s_2':'UC_x_s_1-tau_2','R_y_s_1_y_s_2':'R_x_s_1_y_s_2+lambda_payoff','C_y_s_1_y_s_2':'C_y_s_1_x_s_2+mu_payoff',
 'vR':'R_y_s_1_x_s_2-R_y_1_1_x_s_2','vC':'C_x_s_1_y_s_2-C_x_s_1_y_1_2',
 'Ractual':'Urow-((1-p)*(1-q)*R_x_s_1_x_s_2+p*(1-q)*R_y_s_1_x_s_2+(1-p)*q*R_x_s_1_y_s_2+p*q*R_y_s_1_y_s_2)',
 'Cactual':'Ucol-((1-p)*(1-q)*C_x_s_1_x_s_2+(1-p)*q*C_x_s_1_y_s_2+p*(1-q)*C_y_s_1_x_s_2+p*q*C_y_s_1_y_s_2)'}
assert set(cert['algebraic_definitions'])==set(expected_definitions)
for name,expr in cert['algebraic_definitions'].items():
 equal(parse(expr),parse(expected_definitions[name]))
 definitions[symbols[name]]=parse(expr)
atoms={k:parse(v['expression']) for k,v in cert['nonnegative_atoms'].items()}
premises={k:parse(v['expression']) for k,v in cert['premises'].items()}
for item in list(cert['nonnegative_atoms'].values())+list(cert['premises'].values()):
 assert set(item['provenance'])<=valid_ids
# The initial profile is no longer selected; its lower bound uses the
# implemented OM block's exact input-profile domination, not B01 alone.
assert cert['premises']['bad_g']['provenance']==['T01','O02','B01']

# Independently constrain the range atoms rather than trusting a JSON
# statement that an arbitrary factor is nonnegative.
expected_atoms={'rho':'rho','sigma':'sigma','eta':'eta','eta_bar':'1-eta',
 'xi':'xi','xi_bar':'1-xi','alpha':'alpha','alpha_bar':'1-alpha',
 'gamma':'gamma','gamma_bar':'1-gamma','delta_f':'delta_f',
 'tau_1':'tau_1','tau_2':'tau_2','delta_f_minus_tau_1':'delta_f-tau_1',
 'delta_f_minus_tau_2':'delta_f-tau_2','p':'p','p_bar':'1-p',
 'q':'q','q_bar':'1-q','theta':'theta','theta_bar':'1-theta',
 'mass':'mass','mass_bar':'1-mass'}
assert set(atoms)==set(expected_atoms)
for key,expr in expected_atoms.items():equal(atoms[key],parse(expr))

derived_premises={k:v for k,v in cert['premises'].items() if 'derived_by' in v}
assert set(derived_premises)=={'old_row_cap','old_col_cap'}
available={k:v for k,v in premises.items() if k not in derived_premises}
for key,expr in atoms.items():available['atom:'+key]=expr
checked=[]
guard_checks=[]
for entry in cert['nonnegative_combinations']:
 assert entry['id'] not in available
 total=s.Integer(0)
 for term in entry['combination']:
  if term['slack'] in derived_premises and term['slack'] not in available:
   # Exhaustive real-order cases implement Min(1,Vstar); the input is the
   # already-proved *unclipped* cancellation, plus normalized regret<=1.
   key=term['slack'];rule=derived_premises[key]['derived_by']
   assert rule['rule']=='minimum_after_cancellation'
   side='row' if key=='old_row_cap' else 'col'
   mass='eta' if side=='row' else 'xi';loss='tau_1' if side=='row' else 'tau_2'
   raw='vR' if side=='row' else 'vC';cap='Vr' if side=='row' else 'Vc'
   assert rule==dict(rule='minimum_after_cancellation',cleared='old_'+side+'_unclipped',
    normalized='old_'+side+'_regret_upper',positive_mass=mass+'_lower',mass=mass,cap=cap)
   assert all(rule[x] in available for x in ['cleared','normalized','positive_mass'])
   f=parse(raw+'+'+loss);ms=parse(mass);vs=parse('(1-'+mass+'-r)/'+mass)
   equal(ms*(vs-f),available[rule['cleared']])
   equal(1-f,available[rule['normalized']])
   equal(ms-parse('r'),available[rule['positive_mass']])
   expected_min=s.Min(1,vs)
   equal(parse(cert['delta_f_free_ideal_aliases'][cap]),expected_min)
   # Branch Vstar<=1: cap-f=(cleared slack)/mass>=0 using mass>=r>0.
   # Branch Vstar>=1: cap-f=1-f>=0 using normalized true-regret upper.
   guard_checks.append(dict(output=key,cap_binding=str(expected_min),
    cases=['Vstar<=1: cleared/mass; mass>=r>0','Vstar>=1: normalized upper']))
   available[key]=premises[key]
  assert term['slack'] in available
  coeff=Fraction(term['coefficient']);assert coeff>=0
  multiplier=s.Rational(coeff.numerator,coeff.denominator)
  for factor in term['factors']:
   assert factor in atoms
   multiplier*=atoms[factor]
  total+=multiplier*available[term['slack']]
 target=parse(entry['target']);equal(target,total)
 available[entry['id']]=target
 checked.append(entry['id'])

# Independently verify the key targets, protecting the theorem semantics
# from a certificate which only replaces a target by a different identity.
expected_targets={
 'eta_upper':'1-r-eta','xi_upper':'1-r-xi',
 'old_row_unclipped':'1-eta-r-eta*(vR+tau_1)',
 'old_col_unclipped':'1-xi-r-xi*(vC+tau_2)',
 'row_bridge':'rho*R_y_s_1_y_1_2-G_s+sigma*Vc','col_bridge':'sigma*C_y_1_1_y_s_2-G_s+rho*Vr',
 'common_column_deficit_plane':'eta*sigma*AC+rho*(1-eta)*DC-rho*r+eta*rho-eta*G_s',
 'common_row_deficit_plane':'xi*rho*AR+sigma*(1-xi)*DR-sigma*r+xi*sigma-xi*G_s',
 'generic_payoff_lower':'mass*(new-base)-r+(1-mass)*cap',
 'rectangle_failure':'Fpq-r-delta_f',
 'rectangle_row_budget':'delta_f*(1-q)-((1-p)*(1-q)*delta_f+p*(1-q)*tau_1)',
 'rectangle_col_budget':'delta_f*(1-p)-((1-p)*(1-q)*delta_f+(1-p)*q*tau_2)',
 'row_family_error':'MI+delta_f-Ractual','col_family_error':'MI+delta_f-Cactual',
}
for key,expr in expected_targets.items():equal(available[key],parse(expr))
for expr in cert['zero_identities'].values():equal(parse(expr),s.Integer(0))

# Guarded clipping is finite arithmetic, not a claim that a raw difference
# must itself be nonnegative.  The first rule is applied AFTER cancellation.
eta,xi,r,rho,sigma,G_s,vR,vC,tau_1,tau_2=[symbols[z] for z in
 ('eta','xi','r','rho','sigma','G_s','vR','vC','tau_1','tau_2')]
Vstar=(1-eta-r)/eta;Wstar=(1-xi-r)/xi
equal(eta*(Vstar-parse('vR')-tau_1),available['old_row_unclipped'])
equal(xi*(Wstar-parse('vC')-tau_2),available['old_col_unclipped'])
equal(premises['old_row_cap'],parse('Vr-vR-tau_1'))
equal(premises['old_col_cap'],parse('Vc-vC-tau_2'))
clip_cases=[
 {'guard':'Vstar<=1','cap':'Vstar','positive_divisor':'eta>=r>0','slack':'old_row_unclipped/eta'},
 {'guard':'Vstar>=1','cap':'1','positive_divisor':None,'slack':'1-(vR+tau_1), normalized true regret'},
 {'guard':'Wstar<=1','cap':'Wstar','positive_divisor':'xi>=r>0','slack':'old_col_unclipped/xi'},
 {'guard':'Wstar>=1','cap':'1','positive_divisor':None,'slack':'1-(vC+tau_2), normalized true regret'},
]
# For max(0, bridge), the zero branch uses normalized payoff>=0;
# the positive branch uses the exact cleared bridge below.
equal(rho*(symbols['R_y_s_1_y_1_2']-(G_s-sigma*symbols['Vc'])/rho),available['row_bridge'])
equal(sigma*(symbols['C_y_1_1_y_s_2']-(G_s-rho*symbols['Vr'])/sigma),available['col_bridge'])
assert {x['id'] for x in cert['guarded_rules']}=={
 'minimum_after_cancellation','nonnegative_bridge_clip','family_max_transfer'}

# No accumulating error enters the first/second source payoff template.
# Verify all four actual instantiations of its complementary cap formula.
aliases={k:parse(v) for k,v in cert['delta_f_free_ideal_aliases'].items()}
instantiations=cert['generic_source_instantiations']
assert len(instantiations)==4
expected_source={
 's_1':('eta','lambda_payoff','1-mu_payoff','T_y_s_1'),
 's_2':('xi','mu_payoff','1-lambda_payoff','T_y_s_2'),
 'S':('eta','kappa_1','s_2','T_y_s_1'),
 'T':('xi','kappa_2','s_1','T_y_s_2')}
for item in instantiations:
 mass,base_alias,cap_alias,tuple_id=expected_source[item['alias']]
 assert item['mass']==mass and item['tuple']==tuple_id
 expression=1-parse(base_alias)-(r-(1-parse(mass))*parse(cap_alias))/parse(mass)
 actual_alias=aliases[item['alias']].xreplace({symbols['h_1']:aliases['h_1'],symbols['h_2']:aliases['h_2']})
 equal(actual_alias,s.Min(1,expression))
assert cert['second_envelope_provenance']['same_masses']==['alpha','gamma']
assert cert['second_envelope_provenance']['nearBR_error_variables']==[]
for name,expr in aliases.items():
 forbidden={'delta_f','tau_1','tau_2','tau','g'}
 assert not ({str(z) for z in expr.free_symbols}&forbidden)
# Recursively replace finite aliases: the conclusion is still delta_f-free.
expanded={}
for name,expr in aliases.items():
 expanded[name]=expr.xreplace({symbols[k]:v for k,v in expanded.items() if k in symbols})
 assert not ({str(z) for z in expanded[name].free_symbols}&{'delta_f','tau_1','tau_2','tau','g'})

# Positive common-envelope denominators use mass range, never a nearBR
# error allowance.  In particular xi*row_D >= xi*r>0, and mirror.
alpha,gamma=[symbols[z] for z in ('alpha','gamma')]
equal(xi*aliases['row_D']-xi*r,
 xi*available['gamma_upper']+gamma*(1-xi))
equal(eta*aliases['col_D']-eta*r,
 eta*available['alpha_upper']+alpha*(1-eta))

report=dict(success=True,certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
 imports_producer=False,nonnegative_combinations_verified=len(checked),
 exact_zero_identities_verified=len(cert['zero_identities']),
 nonnegative_factor_atom_count=len(atoms),
 single_joint_G_shift_verified=True,old_regret_cancellation_before_min_clip_verified=True,
 first_mass_bounds_without_dividing_by_one_minus_loss_verified=True,
 cross_bridge_one_joint_budget_verified=True,
 common_deficit_planes_nonnegative_combination_verified=True,
 generic_source_payoff_lower_template_verified=True,
 source_alias_instantiations_verified=4,
 whole_parameter_family_error_budgets_verified=True,
 OM_tau_subtracted_once_before_arc_instantiation_verified=True,
 same_alpha_gamma_envelope_inputs_delta_free=True,
 positive_envelope_denominators_verified=True,
 derived_min_cap_premises_verified=len(guard_checks),
 guarded_clipping_cases=clip_cases,
 scope='Finite arithmetic templates given typed contract instantiations and normalized payoff facts. Full ideal cap/envelope provenance is a separate ledger. No TS/LP implementation or quantifier theorem is proved here.',
 remaining_mathematical_interfaces=cert['trust_boundary'])
(base/'delta_f-lift-ledger-replay.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ['guarded_clipping_cases','remaining_mathematical_interfaces']},ensure_ascii=False))
