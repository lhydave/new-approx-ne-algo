"""Independent exact derivation and Fraction Bernstein replay of the two signs."""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json
import sympy as s

base=Path(__file__).resolve().parent
cert=json.loads((base/'postclip-orientation-certificate.json').read_text())
r,rho=s.symbols('r rho');R,W,Z=s.symbols('R W Z');sigma=1-rho
r0=s.Rational(cert['r_lower'])
chi_1=(rho*(1-r)-r)/rho;chi_2=(sigma*(1-r)-r)/sigma
v0=r*(1+sigma)/(sigma*(1-r))
hc_one=1-r*(1+rho)/sigma+rho*(1-r)/sigma
rawT_one=hc_one-r+v0*(chi_1+chi_2-r)
gap=s.cancel((r-rawT_one)*rho*sigma*sigma*(1-r))
uH=(r-1+r*(1+rho)/sigma)/(rho*(1-r)/sigma)
minus_cnv=s.cancel(-(sigma*(rho-2*r)+rho*sigma*(chi_2-r)*uH)*(1-r))
exprs={'positive_r_minus_rawT_one_at_v0':gap,
       'negative_CNv_at_h_equals_r':minus_cnv}

def bernstein_sequential(power,degrees):
    vals=dict(power)
    for axis,n in enumerate(degrees):
        other=[a for a in range(3) if a!=axis]; new={}
        for fixed in product(*(range(degrees[a]+1) for a in other)):
            for j in range(n+1):
                key=[0,0,0]
                for a,t in zip(other,fixed):key[a]=t
                total=F(0)
                for i in range(j+1):
                    key[axis]=i
                    total+=vals.get(tuple(key),F(0))*F(comb(j,i),comb(n,i))
                key[axis]=j;new[tuple(key)]=total
        vals=new
    return vals

results=[]
for ent in cert['polynomials']:
    po=s.Poly(s.expand(exprs[ent['name']].subs(rho,r+(1-2*r)*W).subs(r,r0+(s.Rational(1,3)-r0)*R)),R,W,Z)
    power={idx:F(cv) for idx,cv in po.terms()}
    assert power=={tuple(idx):F(cv) for idx,cv in ent['power_coefficients']}
    assert list(po.degree_list())==ent['degree']
    bs=bernstein_sequential(power,ent['degree'])
    assert len(bs)==ent['count']
    bad=[idx for idx,cv in bs.items() if cv<=0]
    assert sorted(bad)==sorted(tuple(idx) for idx in ent['nonpositive_indices'])
    if bad:
        assert bad==[(0,1,0)] and ent['degree'][1]==3
        b0,b1,b2=[bs[(0,j,0)] for j in range(3)]
        assert b0>0 and b2>0 and 12*b0*b2-9*b1*b1>0
        expected=ent['grouped_quadratics'][0]
        assert [b0,b1,b2]==[F(expected[k]) for k in ('b0','b1','b2')]
        assert 12*b0*b2-9*b1*b1==F(expected['positive_discriminant_margin'])
        # Group B^3_0,B^3_1,B^3_2 into (1-W) times a positive quadratic.
        # The last W Bernstein coefficient and every remaining r slice are positive.
        assert all(cv>0 for idx,cv in bs.items() if idx[0]>0 or idx[1]==3)
    else:
        assert min(bs.values())>0
    results.append({'name':ent['name'],'coefficient_count':len(bs),'strictly_positive_polynomial':True})
out={'success':True,'subdivisions':0,'independent_symbolic_derivation':True,
     'independent_fraction_bernstein_replay':True,'results':results}
(base/'postclip-orientation-certificate-replay.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
