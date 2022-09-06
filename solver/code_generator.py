#! /usr/bin/env python

"""code_generator.py

Usage:
    code_generator.py FILE...

Options:
    -h, --help    show help
    FILE          json file to split or join

Output is include/butcher_methods.hpp
"""

import json
import os, sys
import sympy as sp

output_file = "include/solver/butcher_methods.hpp"

def STLizer(expr):
  math_to_stl = [(f,sp.Function("std::"+str(f),nargs=1)) for f in (sp.sin,sp.cos,sp.exp)]
  math_to_stl.append( (sp.sqrt,sp.Function("std::sqrt",nargs=1)) )

  for old,new in math_to_stl:
    expr = expr.subs(old,new)
  
  std_pow = sp.Function("std::pow",nargs=2)
  expr = expr.replace(sp.Pow,std_pow)

  return str(expr)

def is_explicit(filename):
  with open(filename,'r') as f:
    data = json.load(f)

  mA = sp.Matrix(data['A'])

  nstages = len(data['A'])
  return ( sum([
      sum([ sp.Abs(mA[i,j]) for j in range(i,nstages) ])
      for i in range(0,nstages)
  ]) == 0 )

def extract_meth(filename):
  with open(filename,'r') as f:
    data = json.load(f)

  mA = sp.Matrix(data['A'])

  nstages = len(data['A'])
  is_explicit = ( sum([
      sum([ sp.Abs(mA[i,j]) for j in range(i,nstages) ])
      for i in range(0,nstages)
  ]) == 0 )

  label = os.path.basename(filename).replace(".json","")
  A = [ [ STLizer(sp.parse_expr(aij)) for aij in ai ] for ai in data['A'] ]
  b = [ STLizer(sp.parse_expr(bi)) for bi in data['b'] ]
  c = [ STLizer(sp.parse_expr(ci)) for ci in data['c'] ]
  
  return {
    'label':label,
    'A':[ " , ".join(ai) for ai in A ],
    'b':" , ".join(b),
    'c':" , ".join(c)
    }
  
if __name__ == '__main__':
    from docopt import docopt

    args = docopt( __doc__ , sys.argv[1:] )

    list_meth = [ extract_meth(f) for f in args['FILE'] if is_explicit(f) ]

    import jinja2
    env = jinja2.Environment(loader=jinja2.FileSystemLoader("."))
    template = env.get_template("tpl_butcher_methods.hxx")

    with open(output_file,'w') as f:
      f.write(
        template.render(list_meth=list_meth)
      )
