#! /usr/bin/env python

"""code_generator.py

Usage:
    code_generator.py FILE... (--output=<header>)

Options:
    -h, --help                      show help
    FILE                            json file to split or join and name of output (last)
    -o=<header>, --output=<header>  name of output file header C++

Output is include/butcher_methods.hpp
"""

import json
import os, sys
import sympy as sp

def code_pow(x,e):
  """
    replace some trivial cases by a simpler expression (for code)
    because all division in sympy are `sp.Pow(expr,-1)` and same
    for every square, of square root
  """
  if e == 1 :
    return x
  if sp.S(e).is_number and e < 0 :
    return 1./(code_pow(x,-e))
  if e == sp.S.Half :
    return sp.Function("std::sqrt",nargs=1)(x)
  return sp.Function("std::pow",nargs=2)(x,e)

def STLizer(expr):
  math_to_stl = [(f,sp.Function("std::"+str(f),nargs=1)) for f in (sp.sin,sp.cos,sp.exp)]
  math_to_stl.append( (sp.sqrt,sp.Function("std::sqrt",nargs=1)) )

  #expr = 1.0*expr
  expr = 1.0*expr

  for old,new in math_to_stl:
    expr = expr.subs(old,new)
  
  expr = expr.replace(sp.Pow,code_pow)

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

  id = os.path.basename(filename).replace(".json","")
  print(id)
  A = [ [ STLizer(sp.parse_expr(aij)) for aij in ai ] for ai in data['A'] ]
  b = [ STLizer(sp.parse_expr(bi)) for bi in data['b'] ]
  c = [ STLizer(sp.parse_expr(ci)) for ci in data['c'] ]

  butcher = {
    'label': data['label'],
    'id': id,
    'A':[ " , ".join(ai) for ai in A ],
    'b':" , ".join(b),
    'c':" , ".join(c)
    }

  if 'b2' in data:
    b2 = [ STLizer(sp.parse_expr(bi)) for bi in data['b2'] ]
    butcher['b2'] = " , ".join(b2)
  
  return butcher

  
  
if __name__ == '__main__':
    from docopt import docopt

    args = docopt( __doc__ , sys.argv[1:] )

    list_meth = [ extract_meth(f) for f in args['FILE'] if is_explicit(f) ]

    local_dir = os.path.dirname(os.path.abspath(__file__))

    import jinja2
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(local_dir))
    template = env.get_template("tpl_butcher_methods.hxx")

    with open(args['--output'],'w') as f:
      f.write(
        template.render(list_meth=list_meth)
      )
