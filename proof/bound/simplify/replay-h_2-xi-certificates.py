"""Independent exact Bernstein reconstruction; does not import the generator."""
import json
from math import comb
from pathlib import Path
import sympy as s

base=Path(__file__).resolve().parent
data=json.loads((base/'h_2-xi-small-certificates.json').read_text())
rho,r,z,w,t=s.symbols('rho r z w t')
symbols={'rho':rho,'r':r,'z':z,'w':w,'t':t}
results=[]
for chart in data['charts']:
    vars=[symbols[k] for k in chart['variables']]
    ns=chart['degree']
    co={tuple(map(int,k.split(','))):s.Rational(v) for k,v in chart['coefficients'].items()}
    assert len(co)==chart['coefficient_count']
    assert all(v>0 for v in co.values())
    assert min(co.values())==s.Rational(chart['minimum'])
    rebuilt=sum(b*s.prod(comb(n,i)*a**i*(1-a)**(n-i) for a,n,i in zip(vars,ns,ids))
                for ids,b in co.items())
    plo,phi=map(s.Rational,chart['rho_range']);rlo,rhi=map(s.Rational,chart['r_range'])
    original=s.sympify(chart['polynomial'],locals=symbols)
    target=original.subs({rho:plo+(phi-plo)*z,r:rlo+(rhi-rlo)*w})
    assert s.Poly(s.expand(rebuilt-target),*vars).is_zero
    results.append({'name':chart['name'],'coefficient_count':len(co),'exact_reconstruction':True})
assert data['symmetric_cap_quartic_at_lower']==str((5*r**4+2*r**2-4*r+1).subs(r,s.Rational(619,2000)))
report={'all_passed':True,'scope':'Independent JSON polynomial reconstruction and positive exact coefficients','charts':results}
(base/'h_2-xi-certificates-replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'all_passed':True,'charts':len(results),'coefficients':sum(k['coefficient_count'] for k in results)}))
