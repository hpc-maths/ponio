#! /usr/bin/env python

"""analysis.py

Usage:
  analysis.py (--input=<input_json>) (--output=<output_dir>) [--verbose]

Options:
  -h, --help                                 Show help
  -i <input_json>, --input=<input_json>      Input file which contains Butcher tableaus
  -o <output_dir>, --output=<output_dir>     Output directory where store Runge-Kutta scheme analysis and static API on stability domain and order stars
  --verbose                                  Verbose mode

"""

import sympy as sp
import numpy as np

import json
import os, sys

class rk_butcher:
    def __init__(self,label,A,b,c):
        self.label = label
        self.A = sp.Matrix(A)
        self.b = sp.Matrix(b)
        self.c = sp.Matrix(c)

        self._order = None
        self._stage_order = None
        self._is_explicit = None
        self._is_approched_method = None
        self._R = None
        self._A = None
        self._P = None
        self._id = None
    
    def __iter__(self):
        """ to convert object into dictionary
        """
        yield 'label',self.label
        yield 'A', [
            [ str(aij) for aij in self.A.row(i) ]
            for i in range(self.nstages)
        ]
        yield 'b', [ str(bi) for bi in self.b ]
        yield 'c', [ str(ci) for ci in self.c ]

    @property
    def nstages(self):
        return len(self.b)
    
    @property
    def order(self):
        r""" compute order of a Runge-Kutta method

        turn around (0,0) (with a radius of `r=0.25`) and count changement of sign of order star function
        """
        if self._order is None:
            self._order = 0
            z = sp.symbols("z")
            r = 0.25
            Z = r*np.exp(1j*np.linspace(0.,2.*np.pi))
            tmp = sp.lambdify(z, sp.Abs(self.order_star_function())-1 , 'numpy' )(Z)
            for i,a in enumerate(tmp[:-1]):
                if a*tmp[i+1] < 0:
                    self._order += 1
            self._order = int(self._order/2)-1
        return self._order
    
    @property
    def stage_order(self):
        r""" compute stage order of a Runge-Kutta method with the following equation:

        The stage order $r$ is the highest number such:
        $$
            \forall i,\ \sum_j a_{ij}c_j^{s-1} = \frac{1}{s}c_i^s,\quad s=1,\dots,r
        $$
        """
        if self._stage_order is None:
            self._stage_order = sum([
                sum([
                    sum([ self.A[i,j]*self.c[j]**(s-1) for j in range(self.nstages) ]) == sp.Rational(1,s)*self.c[i]**s
                    for i in range(self.nstages)
                ]) == self.nstages # test is all equality are true
                for s in range(1,self.order+1) # could be a while...
            ])
        return self._stage_order
    
    @property
    def is_explicit(self):
        if self._is_explicit is None:
            self._is_explicit = ( sum([
                sum([ sp.Abs(self.A[i,j]) for j in range(i,self.nstages) ])
                for i in range(0,self.nstages)
            ]) == 0 )
        return self._is_explicit

    @property
    def is_approched_method(self):
        if self._is_approched_method is None:
            self._is_approched_method = any([ type(x) is sp.Float for x in self.A ])
        return self._is_approched_method

    def stability_function(self):
        if self._R is None:
            z = sp.symbols("z")
            if self.is_approched_method:
                # solve Dahlquist test equation
                # because SymPy can't simplify determinant formula
                """
                stages = self.scheme(latex=False)
                dtsubs = (sp.symbols("dt"),z)
                unsubs = (sp.symbols("un"),1)
                fsubs  = (sp.Function("f",nargs=2),lambda t,u:u)
                self._R = stages[-1].rhs
                
                if not self.is_explicit:
                    Ki = []
                    for eq in stages[:-1]:
                        tmp = eq.replace(*fsubs).subs(*unsubs).subs(*dtsubs)
                        Ki.append( (eq.lhs, sp.solve(tmp,eq.lhs)[0].simplify()) )
                    for ki,sub in reversed(Ki):
                        self._R = self._R.subs(ki,sub)
                    self._R = self._R.subs(*unsubs).subs(*dtsubs).expand()
                else :
                    for eq in reversed(stages[:-1]):
                        self._R = self._R.subs(eq.lhs,eq.rhs)
                    self._R = self._R.replace(*fsubs).subs(*unsubs).subs(*dtsubs).factor().evalf()
                """
                self._R = ( 1+z*( self.b.T*( (sp.eye(*self.A.shape) - z*self.A).inv().evalf(chop=True)*sp.ones(*self.b.shape) ) )[0] ).expand().simplify().collect(z)
            else:
                # standard formula from Butcher's book
                I  = sp.eye(*self.A.shape)  # identity matrix
                I1 = sp.ones(*self.b.shape) # vector of ones
                N = I + z*( I1*self.b.T - self.A )
                D = I - z*self.A
                self._R = ( N.det()/D.det() ).expand().simplify().collect(z)
        return self._R
    
    def order_star_function(self):
        if self._A is None:
            z = sp.symbols("z")
            self._A = sp.exp(-z)*self.stability_function()
        return self._A
    
    def precision_function(self):
        if self._P is None:
            z = sp.symbols("z")
            self._P = sp.Abs( (self.stability_function() - sp.exp(z))/sp.exp(z) )
        return self._P
    
    def scheme(self,latex=True,lawson=False):
        # define a zero
        zero = 0 if not lawson else sp.zeros(3,1)

        if latex:
            ks = sp.symbols("k_{{1:{}}}".format(self.nstages+1))
            un,unp1 = sp.symbols(r"u^n u^{n+1}")
            tn,dt = sp.symbols(r"t^n \Delta\ t")

            # to keep exponential in first, use matrix expressions
            if lawson:
                ks = [ sp.MatrixSymbol(sp.latex(ki),3,1) for ki in ks ]
                un,unp1 = ( sp.MatrixSymbol(sp.latex(ui),3,1) for ui in (un,unp1) )
        else:
            ks = sp.symbols("k1:{}".format(self.nstages+1))
            un,unp1 = sp.symbols(r"un unp1")
            tn,dt = sp.symbols(r"tn dt")

            # to keep exponential in first, use matrix expressions
            if lawson:
                ks = [ sp.MatrixSymbol(str(ki),3,1) for ki in ks ]
                un,unp1 = ( sp.MatrixSymbol(str(ui),3,1) for ui in (un,unp1) )
        
        f = sp.Function("f",nargs=2)
        a = 1
        if lawson:
            def mexp(expr):
                if type(expr) is sp.ZeroMatrix:
                    return 1
                if latex:
                    expr_str = r"e^{{ {} }}".format(sp.latex(expr.simplify()))
                else:
                    expr_str = r"exp({})".format(str(expr.simplify().evalf()).replace("1.0*",""))
                return sp.MatrixSymbol( expr_str , *expr.shape)
            
            L = sp.MatrixSymbol("L",3,3)
            #N = sp.Function("N",nargs=2)
            def N(t,u):
                if latex:
                    expr_str = r"N\left({}, {}\right)".format(sp.latex(t),sp.latex(u.simplify()))
                else:
                    expr_str = r"N({}, {})".format(str(t),str(u.simplify()))
                return sp.MatrixSymbol( expr_str , *u.shape )
            f = lambda t,u:mexp((tn-t)*L)*N(t,mexp((t-tn)*L)*u)
            a = mexp(dt*L) # make the last changement of variable with u_{n+1}

        stages = []
        for i in range(0,self.nstages):
            args = [ self.A[i,j]*ks[j] for j in range(0,self.nstages)]
            if not lawson:
                expr = sum(args)
            else:
                expr = sp.MatAdd(*args)
            stages.append(sp.Eq(
                ks[i],
                f( tn+self.c[i]*dt , un + dt*expr )
            ))
        
        args = [ self.b[i]*ks[i] for i in range(0,self.nstages) ]
        if not lawson:
            expr = sum(args)
        else:
            expr = sp.MatAdd(*args)
        stages.append(sp.Eq(
            unp1,
            a*( un + dt*expr )
        ))
        return stages

    @property
    def id(self):
        if self._id is None:
            self._id = self.label.lower()
            replacements = [
                        (" ","_"),
                        ("(",""),
                        (")",""),
                        (",",""),
                        ("-",""),
                        ("/",""),
                    ]
            for old,new in replacements :
                self._id = self._id.replace(old,new)
        return self._id

    def code(self,lawson=False):
        def func_name(rk_id,lawson=False):
            return ("lawson_" if lawson else "")+rk_id
        
        fun = "f"
        if lawson:
            fun = "L, N"
        args = [fun,"tn","un","dt"] + (["exp=np.exp"] if lawson else [] )
        r = [ "def {funcname}( {args} ):".format( funcname=func_name(self.id,lawson), args=", ".join(args) ) ]
        r.extend([ "{} = {}".format(str(eq.lhs),str(eq.rhs)) for eq in self.scheme(latex=False,lawson=lawson) ])
        r[-1] = r[-1].replace("unp1 =","return")

        return "\n\t".join(r)

def evalf_Cdomain( z , expr , zmin , zmax , N ):
    func = sp.lambdify(z,expr,'numpy')
    X,Y = np.meshgrid(
                np.linspace(np.real(zmin),np.real(zmax),N[0]),
                np.linspace(np.imag(zmin),np.imag(zmax),N[1])
            )
    Z = X + 1j*Y
    return np.abs( func(Z) )

from docopt import docopt
import re

if __name__ == '__main__':
    args = docopt( __doc__ , sys.argv[1:] )
    
    vprint = print if args['--verbose'] else lambda *args,**kwargs:None

    with open( args['--input'] ,'r') as f:
        data = json.load(f)

    output = args['--output']
    os.makedirs(output,exist_ok=True)

    # to remove parenthesis around fractions
    pattern = re.compile("\\\\left\(\\\\frac\{(?P<numerator>[^{}]+)\}\{(?P<denominator>[^{}]+)\}\\\\right\)")

    rk_analysis = []
    for i,drk in enumerate(data):
        vprint("{}/{} {:25}".format(i+1,len(data),drk['label']), end="\n")

        rk = rk_butcher(drk['label'],drk['A'],drk['b'],drk['c'])
        drk['id'] = rk.id
        drk['nstages'] = rk.nstages
        drk['order']   = rk.order
        drk['stage_order']   = rk.stage_order
        drk['stability_function'] = str(rk.stability_function())
        drk['scheme'] = [ "{} = {}".format(sp.latex(eq.lhs),sp.latex(eq.rhs)) for eq in rk.scheme() ]
        drk['lawson_scheme'] = [ "{} = {}".format(sp.latex(eq.lhs),pattern.sub( lambda m:m.group(0)[6:-7], sp.latex(eq.rhs))) for eq in rk.scheme(lawson=True) ]
        if rk.is_explicit:
            drk['code'] = rk.code()
            drk['lawson_code'] = rk.code(lawson=True)
        rk_analysis.append(drk)

        zmin = -5-3.5j
        zmax =  2+3.5j
        N = (512,512)
        stability_domain = { 'xmin':np.real(zmin), 'xmax':np.real(zmax), 'ymin':np.imag(zmin), 'ymax':np.imag(zmax) }
        stability_domain['data'] = evalf_Cdomain( sp.symbols("z") , rk.stability_function() , zmin , zmax , N ).tolist()

        zmin = -3.5-3.5j
        zmax =  3.5+3.5j
        N = (512,512)
        order_star = { 'xmin':np.real(zmin), 'xmax':np.real(zmax), 'ymin':np.imag(zmin), 'ymax':np.imag(zmax) }
        order_star['data'] = evalf_Cdomain( sp.symbols("z") , rk.order_star_function() , zmin , zmax , N ).tolist()

        rk_domain = os.path.join(output,rk.id+".json")
        json.dump( {'stability_domain':stability_domain, 'order_star':order_star} , open(rk_domain,'w') )

    vprint()
    
    json_analysis = os.path.join(output,"analysis.json")
    json.dump(rk_analysis, open(json_analysis,'w')  ,indent=4)

