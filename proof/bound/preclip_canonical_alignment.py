"""Canonical preclip S >= r may force Ck < r: fixed chart test."""
from math import comb
from pathlib import Path
import json
import sympy as s
r,rho,u,v=s.symbols('r rho u v');sigma=1-rho
chi_1=1-r-r/rho;chi_2=1-r-r/sigma;c=1-2*r
hr0=1-r*(1+sigma)/rho;hr1=sigma*(1-r)/rho
S0=hr0-r+u*(chi_2-r);S1=hr1+u*chi_1
vs=s.factor((r-S0)/S1)
cn=rho*rho*chi_1+rho*rho*c*u+(sigma*(rho-2*r)+rho*sigma*(chi_2-r)*u)*v
expr=s.factor(-cn.subs(v,vs))
num,den=s.fraction(expr)
sign=s.sign(den.subs({r:s.Rational(31,100),rho:s.Rational(9,20),u:s.Rational(4,5)}))
a,b,w=s.symbols('a b w');r0=s.Rational('0.30953996')
rr=r0+(s.Rational(1,3)-r0)*a; pp=rr+(s.Rational(1,2)-rr)*b
uu=rr/(1-rr)+(1-rr/(1-rr))*w
cleared=s.factor((num*sign).subs({r:rr,rho:pp,u:uu}))
n,d=s.fraction(cleared)
pol=s.Poly(s.expand(n),a,b,w);degrees=pol.degree_list();vals=[]
for i in range(degrees[0]+1):
 for j in range(degrees[1]+1):
  for k in range(degrees[2]+1):
   vals.append(sum(co*s.Rational(comb(i,h),comb(degrees[0],h))
     *s.Rational(comb(j,l),comb(degrees[1],l))
     *s.Rational(comb(k,t),comb(degrees[2],t))
     for (h,l,t),co in pol.terms() if h<=i and l<=j and t<=k))
out={'expression':str(expr),'positive_denominator_to_justify':str(s.factor(den*sign)),
     'substitution_denominator':str(d),'degree':list(degrees),'count':len(vals),
     'all_positive':all(z>0 for z in vals),'minimum':str(min(vals)),
     'negative_count':sum(1 for z in vals if z<0)}
print(json.dumps(out,indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
