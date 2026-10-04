"""Independent coefficient replay: direct monomial-to-Bernstein evaluation.
Does not import the producer or use its axis conversion routine.
"""
import sympy as s,json
from pathlib import Path
from math import comb
from itertools import product
from fractions import Fraction as F
path=Path(__file__).with_name('cross_derivative_large_mass_certificate.json');raw=json.loads(path.read_text())
vars=s.symbols('z_r z_rho z_eta z_xi');loc={str(z):z for z in vars};out=[]
for chart in raw['charts']:
 assert chart['certified']
 poly=s.Poly(s.sympify(chart['polynomial'],locals=loc)-s.Rational(chart['multiplier'])*s.sympify(chart['feasibility_polynomial'],locals=loc),*vars)
 co={ind:F(c) for ind,c in poly.terms()};deg=tuple(chart['degree']);calculated=[]
 for j in product(*(range(d+1) for d in deg)):
  val=F()
  for powers,c in co.items():
   if all(k<=a for k,a in zip(powers,j)):
    factor=F(1)
    for a,k,d in zip(j,powers,deg):factor*=F(comb(a,k),comb(d,k))
    val+=c*factor
  assert val>0
  calculated.append(val)
 assert [str(c) for c in calculated]==chart['coefficients']
 assert str(min(calculated))==chart['minimum']
 out.append({'chart':chart['chart'],'coefficient_count':len(calculated),'minimum':str(min(calculated)),'passed':True})
result={'scope':'Independent replay of the free-r two-chart crossing-derivative certificate, not the final approximation-bound proof.','total_coefficients':sum(a['coefficient_count'] for a in out),'charts':out}
Path(__file__).with_name('cross-derivative-large-mass-independent-replay.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
