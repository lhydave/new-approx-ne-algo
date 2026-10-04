"""Proof-object replay for ideal payoff/cap provenance (no game search).

Trusted input: the explicitly named Stage-I K/T/R/C/O/B logical contracts.
Checked work: typed same-game instances, universal instantiation, nonnegative
combinations, selector cases, affine mixing, identities and positive divisions.
"""
from pathlib import Path
import argparse,json,hashlib
import sympy as s

HERE=Path(__file__).resolve().parent

class Replay:
    def __init__(self,data):
        self.data=data;assert data['schema']==1 and data['delta']==0
        assert data['game']=='G' and not data.get('extra_building_block_properties')
        self.rows=set(data['row_strategies']);self.cols=set(data['column_strategies'])
        assert {'x_s_1','y_s_1','y_1_1','y_2_1','h','x_1_1'}<=self.rows
        assert {'x_s_2','y_s_2','y_1_2','y_2_2','k','x_1_2'}<=self.cols
        self.names={z:s.Symbol(z)for z in data['symbols']}
        self.terms=data['payoff_terms']
        for name,term in self.terms.items():
            assert term['game']=='G'
            if term['player']in('R','C'):
                row,col=term['arguments'];assert row in self.rows and col in self.cols
                assert name==term['player']+'_'+row+'_'+col
            elif term['player']=='UR':
                assert len(term['arguments'])==1 and term['arguments'][0]in self.cols
                assert name=='UR_'+term['arguments'][0]
            else:
                assert term['player']=='UC' and len(term['arguments'])==1 and term['arguments'][0]in self.rows
                assert name=='UC_'+term['arguments'][0]
        self.defs={self.names[n]:self.parse(e)for n,e in data['rational_definitions'].items()}
        self.known={};self.kinds={};self.positive=[];self.nonnegative=[]
        for name in ('rho','sigma','r'):self.positive.append(self.parse(name))
        for name in ('eta','xi','alpha','gamma','pMix','qMix'):
            self.nonnegative.extend([self.parse(name),1-self.parse(name)])
        assert self.parse(data['assumptions']['front_mass_identity'])==self.parse('rho+sigma-1')
        assert data['assumptions']['positive']==['rho','sigma','r']
        self.tuples=data['tuples']
        assert self.tuples=={
            'y_s_1':{'clause':'R02','input':'y_s_1','response':'y_1_1','mass':'eta','value':'Hw'},
            'y_s_2':{'clause':'C02','input':'y_s_2','response':'y_1_2','mass':'xi','value':'Hz'},
            'y_1_1':{'clause':'R02','input':'y_1_1','response':'y_2_1','mass':'alpha','value':'Hr'},
            'y_1_2':{'clause':'C02','input':'y_1_2','response':'y_2_2','mass':'gamma','value':'Hj'}}
        self.check_definitions();self.counts={}

    def parse(self,expr):return s.sympify(expr,locals=self.names)
    def expand(self,expr):
        for repeat in range(5):
            out=expr.xreplace(self.defs)
            if out==expr:break
            expr=out
        return s.cancel(expr)
    def equal(self,a,b):return s.cancel(self.expand(a-b))==0
    def R(self,x,y):return self.names['R_'+x+'_'+y]
    def C(self,x,y):return self.names['C_'+x+'_'+y]
    def UR(self,y):return self.names['UR_'+y]
    def UC(self,x):return self.names['UC_'+x]

    def sign(self,expr,strict=False):
        if expr==0:return not strict
        num,den=s.fraction(s.factor(expr));sign=1
        for part,isden in ((num,False),(den,True)):
            coefficient,factors=s.factor_list(part)
            if coefficient==0:return not strict and not isden
            if coefficient<0:sign*=-1
            for factor,power in factors:
                pool=self.positive if isden or strict else self.positive+self.nonnegative
                ok=False
                for base in pool:
                    if s.factor(factor-base)==0:ok=True;break
                    if s.factor(factor+base)==0:
                        if power%2:sign*=-1
                        ok=True;break
                if not ok:return False
        return sign>0

    def check_definitions(self):
        n=self.names;expected={
            'lambda_0':'g/rho','mu_0':'g/sigma','u':'(1-eta)/eta','v':'(1-xi)/xi',
            'h_1':'1-kappa_1','h_2':'1-kappa_2','DR':'1-gamma+gamma*v','DC':'1-alpha+alpha*u',
            'LC':'r-eta*(1-lambda_0)','LR':'r-xi*(1-mu_0)',
            'ER':'((1-gamma)*S+gamma*(1-r/xi+v*s_1)-r)/DR',
            'EC':'((1-gamma)*T+gamma*v*h_2)/DR',
            'CR':'((1-alpha)*S+alpha*u*h_1)/DC',
            'CC':'((1-alpha)*T+alpha*(1-r/eta+u*s_2)-r)/DC',
            'JCexact':'(1-pMix)*UC_x_s_1+pMix*UC_y_s_1-UC_xMix',
            'JRexact':'(1-qMix)*UR_x_s_2+qMix*UR_y_s_2-UR_yMix',
            'Rcap':'(1-pMix)*(1-qMix)*g+qMix*(1-pMix*lambda_0)-kR*LR',
            'Ccap':'(1-pMix)*(1-qMix)*g+pMix*(1-qMix*mu_0)-kC*LC'}
        assert set(self.data['rational_definitions'])==set(expected)
        for key,value in expected.items():assert self.equal(n[key],self.parse(value)),key
        clips=self.data['selectors']
        expected_clips={
            'vr':('min',['1','(1-eta-r)/eta']),
            'vc':('min',['1','(1-xi-r)/xi']),
            's_1':('min',['1','1-lambda_0-(r-(1-eta)*(1-mu_0))/eta']),
            's_2':('min',['1','1-mu_0-(r-(1-xi)*(1-lambda_0))/xi']),
            'kappa_1':('max',['0','lambda_0-sigma*vc/rho']),
            'kappa_2':('max',['0','mu_0-rho*vr/sigma']),
            'S':('min',['1','1-kappa_1-(r-(1-eta)*s_2)/eta']),
            'T':('min',['1','1-kappa_2-(r-(1-xi)*s_1)/xi']),
            'kC':('min',['pMix/(1-eta)','(1-pMix)*rho/(eta*sigma)']),
            'kR':('min',['qMix/(1-xi)','(1-qMix)*sigma/(xi*rho)'])}
        assert set(clips)==set(expected_clips)
        for key,(op,values)in expected_clips.items():
            assert clips[key]['operation']==op
            assert len(clips[key]['choices'])==2
            for x,y in zip(clips[key]['choices'],values):assert self.equal(self.parse(x),self.parse(y))
        mixes=self.data['mixtures']
        assert set(mixes)=={'yE','xE','xMix','yMix'}
        expected_mixes={'yE':('column','y_1_2','y_s_2','gamma*v/DR'),
            'xE':('row','y_1_1','y_s_1','alpha*u/DC'),
            'xMix':('row','x_s_1','y_s_1','pMix'),'yMix':('column','x_s_2','y_s_2','qMix')}
        for name,(side,left,right,beta)in expected_mixes.items():
            mix=mixes[name];assert (mix['type'],mix['left'],mix['right'])==(side,left,right)
            assert self.equal(self.parse(mix['coefficient']),self.parse(beta))
            assert mix['game']=='G'

    def primitive(self,node):
        op=node['operation'];n=self.names;r=n['r']
        if op=='payoff_bound':
            assert node['clause_ids']==['K01']
            v=self.parse(node['quantity']);assert str(v)in self.terms
            assert self.terms[str(v)]['player']in('R','C')
            return v if node['side']=='lower'else 1-v,'ge'
        if op=='maximum_bound':
            assert node['clause_ids']==['K01','K02']
            v=self.parse(node['quantity']);assert str(v)in self.terms and self.terms[str(v)]['player']in('UR','UC')
            return 1-v,'ge'
        if op=='BR_upper':
            assert node['clause_ids']==['K02']
            if node['player']=='R':return self.UR(node['opponent'])-self.R(node['deviation'],node['opponent']),'ge'
            return self.UC(node['opponent'])-self.C(node['opponent'],node['deviation']),'ge'
        if op=='plane':
            key=node['tuple'];binding=self.tuples[key];m=n[binding['mass']];H=n[binding['value']]
            assert node['clause_ids']==[binding['clause']]
            x,y=binding['input'],node['opponent'];response=binding['response']
            if binding['clause']=='R02':
                assert y in self.cols
                return m*(self.R(response,y)-self.R(x,y))+(1-m)*(self.UC(x)-self.C(x,y))-H,'ge'
            assert y in self.rows
            return (1-m)*(self.UR(x)-self.R(y,x))+m*(self.C(y,response)-self.C(y,x))-H,'ge'
        if op=='bad_value':
            # g is not selected directly: O02 gives g >= F(OM), and
            # B01 makes that same OM candidate high in a bad execution.
            assert node['clause_ids']==['B01','O02','T01']if node['value']=='g'else node['clause_ids']==['B01']
            assert node['value']in('g','Hw','Hz','Hr','Hj')
            return n[node['value']]-r,'ge'
        if op=='front_BR':
            assert node['clause_ids']==['T02','K02','K01']
            if node['player']=='R':return self.R('y_s_1','x_s_2')-self.UR('x_s_2'),'eq'
            return self.C('x_s_1','y_s_2')-self.UC('x_s_1'),'eq'
        if op=='equal_regret':
            assert node['clause_ids']==['T01','K03']
            if node['player']=='R':return self.UR('x_s_2')-self.R('x_s_1','x_s_2')-n['g'],'eq'
            return self.UC('x_s_1')-self.C('x_s_1','x_s_2')-n['g'],'eq'
        if op=='attainment':
            assert node['clause_ids']==['K02']
            witness=node['witness'];bind=self.data['witnesses'][witness]
            assert bind['game']=='G' and bind['input']==node['opponent'] and bind['player']==node['player']
            if node['player']=='R':return self.R(witness,node['opponent'])-self.UR(node['opponent']),'eq'
            return self.C(node['opponent'],witness)-self.UC(node['opponent']),'eq'
        if op=='alignment':
            if node['tuple']=='y_s_1':
                assert 'mass_eta_low'in self.known and node['clause_ids']==['R03','K02']
                assert self.equal(self.known['mass_eta_low'],n['eta']-r)
                assert self.sign(n['eta'],strict=True)
                return self.R('y_1_1','x_1_2')-self.UR('x_1_2'),'eq'
            assert node['tuple']=='y_s_2'and'mass_xi_low'in self.known and node['clause_ids']==['C03','K02']
            assert self.equal(self.known['mass_xi_low'],n['xi']-r)
            assert self.sign(n['xi'],strict=True)
            return self.C('x_1_1','y_1_2')-self.UC('x_1_1'),'eq'
        if op=='Tdual_row':
            assert node['clause_ids']==['T04','T02','K01'];x=node['opponent']
            return -n['rho']*self.R(x,'x_s_2')+n['sigma']*(self.C(x,'y_s_2')-self.C(x,'x_s_2'))-n['A_s'],'ge'
        if op=='Tdual_col':
            assert node['clause_ids']==['T04','T02','K01'];y=node['opponent']
            return -n['sigma']*self.C('x_s_1',y)+n['rho']*(self.R('y_s_1',y)-self.R('x_s_1',y))-n['B_s'],'ge'
        if op=='Tstop':assert node['clause_ids']==['T05'];return n['V']-n['g'],'ge'
        if op=='Tvalue':
            assert node['clause_ids']==['T05']
            return n['V']-n['rho']*self.R('x_s_1','x_s_2')-n['sigma']*self.C('x_s_1','x_s_2')-n['A_s']-n['B_s'],'eq'
        raise AssertionError(op)

    def run(self):
        for node in self.data['proof']:
            ident=node['id'];assert ident not in self.known
            expr=self.parse(node['expression']);kind=node['kind'];relation=node['relation']
            if kind=='primitive':
                expected,rel=self.primitive(node);assert relation==rel and self.equal(expr,expected),ident
            elif kind=='linear':
                total=s.Integer(0)
                for term in node['terms']:
                    ref=term['fact'];assert ref in self.known;weight=self.parse(term['weight'])
                    if self.kinds[ref]=='ge':
                        if 'sign_fact'in term:
                            signref=term['sign_fact'];assert signref in self.known and self.kinds[signref]=='ge'
                            assert not(self.known[signref].free_symbols&{self.names[z]for z in self.terms})
                            assert self.equal(weight,self.known[signref])
                        else:assert self.sign(weight),('weight not certified',ident,str(weight))
                    total+=weight*self.known[ref]
                assert self.equal(total,expr),('identity failed',ident,s.factor(self.expand(total-expr)))
                assert relation=='ge'
            elif kind=='positive_alias':
                assert relation=='ge'and self.equal(expr,self.names[node['alias']])
                parts=[self.parse(z)for z in node['parts']]
                assert self.equal(expr,sum(parts))and all(self.sign(z)for z in parts)and any(self.sign(z,True)for z in parts)
                self.positive.append(self.names[node['alias']])
            elif kind=='mixing':
                mix=self.data['mixtures'][node['mixture']];beta=self.parse(mix['coefficient'])
                left=self.parse(mix['left_coefficient'])
                assert self.equal(left,1-beta)and self.sign(beta)and self.sign(left)
                other=node['other'];player=node['player'];fn=self.R if player=='R'else self.C
                if mix['type']=='column':expected=fn(other,node['mixture'])-(1-beta)*fn(other,mix['left'])-beta*fn(other,mix['right'])
                else:expected=fn(node['mixture'],other)-(1-beta)*fn(mix['left'],other)-beta*fn(mix['right'],other)
                assert relation=='eq'and self.equal(expr,expected)
            elif kind=='selector':
                assert relation=='ge';alias=self.names[node['alias']];definition=self.data['selectors'][node['alias']]
                choices=[self.parse(z)for z in definition['choices']]
                refs=node['premises']
                for branch,choice in enumerate(choices):
                    target=expr.subs(alias,choice)
                    guard=(choices[1-branch]-choice)if definition['operation']=='min'else(choice-choices[1-branch])
                    valid=target==0 or self.sign(target)or self.equal(target,guard)
                    valid=valid or any(self.equal(target,self.known[ref])and self.kinds[ref]=='ge'for ref in refs)
                    assert valid,('selector branch',ident,branch,str(target))
            elif kind=='instantiate':
                source=node['fact'];assert source in self.known;variable=node['variable'];replacement=node['replacement']
                assert variable in ('h','k')
                assert replacement in(self.rows if variable=='h'else self.cols)
                substitutions={}
                for name,term in self.terms.items():
                    if variable in term['arguments']:
                        args=[replacement if z==variable else z for z in term['arguments']]
                        target=term['player']+'_'+'_'.join(args)
                        assert target in self.names
                        substitutions[self.names[name]]=self.names[target]
                assert self.equal(expr,self.known[source].xreplace(substitutions))and relation==self.kinds[source]
            else:raise AssertionError(kind)
            self.known[ident]=expr;self.kinds[ident]=relation;self.counts[kind]=self.counts.get(kind,0)+1
            if relation=='ge':
                for name in ('g','eta','1-eta','xi','1-xi','alpha','1-alpha','gamma','1-gamma'):
                    z=self.parse(name)
                    if self.equal(expr,z-self.names['r']):self.positive.append(z)
                if not(expr.free_symbols&{self.names[z]for z in self.terms}):self.nonnegative.append(expr)
        for output in self.data['outputs'].values():assert output in self.known and self.kinds[output]=='ge'
        return self.counts

def main():
    ap=argparse.ArgumentParser();ap.add_argument('ledger',nargs='?',default=str(HERE/'ideal-cap-provenance-ledger.json'));args=ap.parse_args()
    path=Path(args.ledger);assert path.stat().st_size<=500000
    data=json.loads(path.read_text());checker=Replay(data);counts=checker.run()
    result={'success':True,'delta_zero_only':True,'same_game_terms_typechecked':True,
        'primitive_clause_instances_rebuilt':True,'nonnegative_multipliers_checked':True,
        'positive_divisions_checked':True,'all_identities_exact':True,'selector_cases_checked':True,
        'proof_node_count':len(data['proof']),'counts':counts,'outputs':data['outputs'],
        'extra_building_block_properties':False,'g_to_r_projection_proved_here':False,
        'ledger_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out=HERE/'ideal-cap-provenance-ledger-replay.json';out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
