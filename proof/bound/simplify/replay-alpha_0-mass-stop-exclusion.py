"""Independent JSON reconstruction; does not import the generating script."""
import json
from math import comb
from pathlib import Path
import sympy as s

base=Path(__file__).resolve().parent
data=json.loads((base/'alpha_0-mass-stop-exclusion.json').read_text())
rho,r,u,v,z,w,t=s.symbols('rho r u v z w t')
symbols={str(k):k for k in [rho,r,u,v,z,w,t]}
results=[]
for chart in data['charts']:
    variables=[symbols[k] for k in chart['variables']];ns=chart['degree']
    mapping={symbols[k]:s.sympify(value,locals=symbols) for k,value in chart['mapping'].items()}
    co={tuple(map(int,k.split(','))):s.Rational(value) for k,value in chart['coefficients'].items()}
    assert len(co)==chart['coefficient_count']
    assert all(value>0 for value in co.values())
    assert min(co.values())==s.Rational(chart['minimum'])
    rebuilt=sum(b*s.prod(comb(n,i)*a**i*(1-a)**(n-i) for a,n,i in zip(variables,ns,ids)) for ids,b in co.items())
    original=s.sympify(chart['polynomial'],locals=symbols).subs(mapping)
    assert s.Poly(s.expand(rebuilt-original),*variables).is_zero
    results.append({'name':chart['name'],'coefficient_count':len(co),'exact_reconstruction':True})
report={'all_passed':True,'scope':'Independent positive rational coefficient reconstruction','charts':results}
(base/'alpha_0-mass-stop-exclusion-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'all_passed':True,'charts':len(results),'coefficients':sum(k['coefficient_count'] for k in results)}))
