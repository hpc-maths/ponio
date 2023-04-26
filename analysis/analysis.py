# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#! /usr/bin/env python

"""analysis.py

Usage:
  analysis.py (FILE...) (--output=<output_dir>) [--verbose] [--standalone]

Options:
  -h, --help                                 Show help
  FILE                                       File which contain a Butcher tableau
  -o <output_dir>, --output=<output_dir>     Output directory where store Runge-Kutta scheme analysis and static API on stability domain and order stars
  --verbose                                  Verbose mode
  -s, --standalone                           All analysis in one file (by default each analysis of each method is in its own file)

"""

import sympy as sp
import numpy as np
import itertools

import json
import os, sys

class rk_butcher:
    def __init__(self,label,A,b,c,b2=None):
        self.label = label
        self.A  = sp.Matrix(A)
        self.b  = sp.Matrix(b)
        self.c  = sp.Matrix(c)
        self.b2 = None if b2 is None else sp.Matrix(b2)

        self._order = None
        self._stage_order = None
        self._is_explicit = None
        self._is_dirk     = None
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

        if self.b2 is not None:
            yield 'b2', [ str(bi) for bi in self.b2 ]

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
    def is_dirk(self):
        if self._is_dirk is None:
            self._is_dirk = ( sum([
                sum([ sp.Abs(self.A[i,j]) for j in range(i+1,self.nstages) ])
                for i in range(0,self.nstages)
            ]) == 0 ) and ( sum([ sp.Abs(self.A[i,i]) for i in range(0,self.nstages) ]) != 0 )
        return self._is_dirk

    @property
    def is_approched_method(self):
        if self._is_approched_method is None:
            self._is_approched_method = any([ type(x) is sp.Float for x in itertools.chain(self.A,self.b,self.c) ])
        return self._is_approched_method

    @property
    def is_embedded(self):
        return self.b2 is not None

    def stability_function(self,embedded=False):
        def stab_func(A,b,approched_method=False):
            z = sp.symbols("z")
            I  = sp.eye(*A.shape)  # identity matrix
            I1 = sp.ones(*b.shape) # vector of ones
            if approched_method:
                return  ( 1+z*( b.T*( (I - z*A).inv().evalf(chop=True)*I1 ) )[0] ).expand().simplify().collect(z)
            elif A.shape[0] >= 10:
                return ( 1+z*( b.T*( (I - z*A).inv()*I1 ) )[0] ).expand().simplify().collect(z)
            else:
                N  = I + z*( I1*b.T - A )
                D  = I - z*A

                return ( N.det()/D.det() ).expand().simplify().collect(z)
        if embedded:
            return stab_func(self.A,self.b2,approched_method=self.is_approched_method)

        if self._R is None:
            self._R = stab_func(self.A,self.b,approched_method=self.is_approched_method)
        return self._R

    def order_star_function(self):
        if self._A is None:
            z = sp.symbols("z")
            self._A = sp.exp(-z)*self.stability_function()
        return self._A

    def relative_error_function(self):
        if self._P is None:
            z = sp.symbols("z")
            self._P = sp.Abs( (self.stability_function() - sp.exp(z))/sp.exp(z) )
        return self._P

    def scheme(self,latex=True,lawson=False):
        if latex:
            ks = sp.symbols("k_{{1:{}}}".format(self.nstages+1))
            un,unp1,ubnp1 = sp.symbols(r"u^n u^{n+1} \tilde{u}^{n+1}")
            tn,dt = sp.symbols(r"t^n \Delta\ t")

            # to keep exponential in first, use matrix expressions
            if lawson:
                ks = [ sp.MatrixSymbol(sp.latex(ki),3,1) for ki in ks ]
                un,unp1,ubnp1 = ( sp.MatrixSymbol(sp.latex(ui),3,1) for ui in (un,unp1,ubnp1) )
        else:
            ks = sp.symbols("k1:{}".format(self.nstages+1))
            un,unp1,ubnp1 = sp.symbols(r"un unp1 ubnp1")
            tn,dt = sp.symbols(r"tn dt")

            # to keep exponential in first, use matrix expressions
            if lawson:
                ks = [ sp.MatrixSymbol(str(ki),3,1) for ki in ks ]
                un,unp1,ubnp1 = ( sp.MatrixSymbol(str(ui),3,1) for ui in (un,unp1,ubnp1) )

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

        if self.is_embedded:
            args = [ self.b2[i]*ks[i] for i in range(0,self.nstages) ]
            if not lawson:
                expr = sum(args)
            else:
                expr = sp.MatAdd(*args)
            stages.append(sp.Eq(
                ubnp1,
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

        args = {"f": "function f(t,u)"}
        if lawson:
            args = {"L": "linar part, a scalar or a matrix", "N": "non-linear part, a function N(t,u)"}
        args.update(
            [
                ("tn","current time"),
                ("un","current solution at time tn"),
                ("dt", "time step to the next time")
            ] + ([("exp=np.exp","exponential function to compute exponential of linear part")] if lawson else [] )
        )
        if self.is_embedded:
            args['tol=1e-4'] = "tolerance of embedded time step method"

        # stages
        scheme = self.scheme(latex=False,lawson=lawson)

        # function declaration
        r = [ "def {funcname}( {args} ):".format( funcname=func_name(self.id,lawson), args=", ".join(args.keys()) ) ]

        args_len = len(max(args.keys(),key=len))
        # docstring
        docstring = """\"\"\"{meth.label}
    {embedded}{lawson}Runge-Kuta method with {meth.nstages} stages, of order {meth.order}
    {arguments}

    return {output}

    {warning}
    source: https://josselin.massot.gitlab.labos.polytechnique.fr/ponio/viewer.html#{meth.id}
\"\"\"""".format(
            meth=self,
            embedded = "" if not self.is_embedded else "Embedded ",
            lawson = "" if not lawson else "Lawson ",
            arguments="\n    ".join([ f"{{k:{args_len}}}    {{v}}".format(k=k,v=v) for k,v in args.items() ]),
            #output   ="{} solution at time tn".format(str(scheme[-1].lhs)) if not self.is_embedded else "({unp1}, {ubnp1}) tuple of solutions at time tn of order {order} and {orderbis}".format(unp1=str(scheme[-2].lhs),ubnp1=str(scheme[-1].lhs),order=self.order,orderbis=self.order-1),
            output   ="{} solution at time tn".format(str(scheme[-1].lhs)) if not self.is_embedded else "(t,u,dt,error) tuple where t is the time of the {unp1}, {ubnp1}) tuple of solutions at time tn of order {order} and {orderbis}".format(unp1=str(scheme[-2].lhs),ubnp1=str(scheme[-1].lhs),order=self.order,orderbis=self.order-1),
            warning  ="WARNING: np.exp is only for a scalar linear case, otherwise need to use scipy.linalg.expm for real matrix exponential" if lawson else ""
        )
        r.extend( docstring.split("\n") )

        r.extend([ "{} = {}".format(str(eq.lhs),str(eq.rhs)) for eq in scheme ])

        if self.is_embedded:
            error_computation = """
# computation error

tmp_unp1,tmp_ubnp1 = (np.atleast_1d({unp1}), np.atleast_1d({ubnp1}))
error = np.sqrt( 1./len(tmp_unp1)*sum([
    ( (unp1_i - ubnp1_i)/(1.0 + max(un_i,unp1_i)) )**2
    for un_i,unp1_i,ubnp1_i in zip(np.atleast_1d(un).flat, tmp_unp1.flat, tmp_ubnp1.flat)
]) )
ndt = dt*(tol/error)**({power})
if error > tol:
    return (tn,un,ndt,error)
else:
    return (tn+dt,{unp1},ndt,error)""".format(unp1=str(scheme[-2].lhs),ubnp1=str(scheme[-1].lhs),power=1./self.order)
            r.extend(error_computation.split("\n"))
            #r.append("return ({unp1}, {ubnp1})".format(unp1=str(scheme[-2].lhs),ubnp1=str(scheme[-1].lhs)))
        else:
            r[-1] = r[-1].replace("{unp1} =".format(unp1=str(scheme[-1].lhs)),"return")

        return "\n    ".join(r)

def analysis_butcher(rk):
    # to remove parenthesis around fractions in Lawson scheme expressions
    pattern = re.compile("\\\\left\(\\\\frac\{(?P<numerator>[^{}]+)\}\{(?P<denominator>[^{}]+)\}\\\\right\)")

    analysis = {}
    # update butcher with new values
    analysis['id']          = rk.id
    analysis['nstages']     = rk.nstages
    analysis['order']       = rk.order
    analysis['stage_order'] = rk.stage_order

    # flags
    analysis['is_explicit'] = rk.is_explicit
    analysis['is_dirk']     = rk.is_dirk
    analysis['is_embedded'] = rk.is_embedded

    # some LaTeX
    analysis['stability_function'] = str(rk.stability_function())
    analysis['scheme'] = [
        "{} = {}".format(sp.latex(eq.lhs), sp.latex(eq.rhs))
        for eq in rk.scheme()
    ]
    analysis['lawson_scheme'] = [
        "{} = {}".format(sp.latex(eq.lhs), pattern.sub( lambda m:m.group(0)[6:-7], sp.latex(eq.rhs)))
        for eq in rk.scheme(lawson=True)
    ]

    # some code
    if rk.is_explicit:
        analysis['code']        = rk.code()
        analysis['lawson_code'] = rk.code(lawson=True)

    return analysis


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

    output = args['--output']
    os.makedirs(output,exist_ok=True)

    rk_analysis = []
    for i,file in enumerate(args['FILE']):
        with open(file,'r') as f:
            butcher = json.load(f)
        vprint("{}/{} {:25}".format(i+1,len(args['FILE']),butcher['label']), end="\r")

        rk = rk_butcher(
            label=butcher['label'],
            A=butcher['A'],b=butcher['b'],c=butcher['c'],
            b2=butcher['b2'] if 'b2' in butcher else None
        )
        butcher.update(analysis_butcher(rk))

        if args['--standalone']:
            rk_analysis.append(butcher)
        else:
            with open(os.path.join(output,rk.id+".json"),'w') as f:
                json.dump(butcher,f,indent=4)

        # add some data evaluated on complex domain
        # stability domaine
        zmin = -5-3.5j
        zmax =  2+3.5j
        N = (512,512)
        stability_domain = { 'xmin':np.real(zmin), 'xmax':np.real(zmax), 'ymin':np.imag(zmin), 'ymax':np.imag(zmax) }
        stability_domain['data'] = evalf_Cdomain( sp.symbols("z") , rk.stability_function() , zmin , zmax , N ).tolist()

        if rk.is_embedded:
            stability_domain['data_embedded'] = evalf_Cdomain( sp.symbols("z") , rk.stability_function(embedded=True) , zmin , zmax , N ).tolist()

        # relative error (or precision) thanks to Laurent François
        relative_error = { 'xmin':np.real(zmin), 'xmax':np.real(zmax), 'ymin':np.imag(zmin), 'ymax':np.imag(zmax) }
        relative_error['data'] = evalf_Cdomain( sp.symbols("z") , rk.relative_error_function() , zmin , zmax , N ).tolist()

        # order star
        zmin = -3.5-3.5j
        zmax =  3.5+3.5j
        N = (512,512)
        order_star = { 'xmin':np.real(zmin), 'xmax':np.real(zmax), 'ymin':np.imag(zmin), 'ymax':np.imag(zmax) }
        order_star['data'] = evalf_Cdomain( sp.symbols("z") , rk.order_star_function() , zmin , zmax , N ).tolist()

        with open(os.path.join(output,rk.id+"_domains.json"),'w') as f:
            json.dump(
                {
                    'stability_domain':stability_domain,
                    'order_star':order_star,
                    'relative_error':relative_error
                },
                f
            )

    vprint()

    if args['--standalone']:
        json_analysis = os.path.join(output,"analysis.json")
        json.dump(rk_analysis, open(json_analysis,'w')  ,indent=4)
