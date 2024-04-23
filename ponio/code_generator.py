# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#! /usr/bin/env python

"""code_generator.py

usage: code_generator.py [-h] [-o OUTPUT] [--Ndigit NDIGIT] [FILE ...]

code generator of Runge-Kutta method from them Butcher tableau

positional arguments:
  FILE                  Files which contain a Butcher tableau

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Name of output file header C++
  --Ndigit NDIGIT       number of digit in output Butcher tableau [default: 36]
"""

import sympy as sp
import numpy as np
import itertools

import json
import os
import sys
import argparse


def phi(i, j=None, c=None):
    # set True if j and c are not None (so compute phi_{ij})
    ij = j is not None and c is not None

    classname = f"phi_{i}" if not ij else f"phi_{i}_{j}"

    @classmethod
    def phi_i_eval(cls, z):
        pass

    def phi_i_doit(self, deep=False, **hints):
        z = self.cj*self.args[0]
        if i == 0:
            return sp.exp(z)
        elif z.is_zero:
            return sp.Rational(1, sp.factorial(i))
        else:
            return ((phi(i-1)(z) - phi(i-1)(0))/z).simplify().expand().simplify()

    def phi_i_latex(self, printer):
        return f"\\varphi_{{{i}}}"
        # return str( type(self).mro()[0] ).replace(" ","\ ")

    def phi_ij_latex(self, printer):
        return f"\\varphi_{{{i},{j}}}"

    return type(
        classname,
        (sp.Function,),
        {
            'cj': 1 if not ij else c[j-1],
            'eval':   phi_i_eval,
            'doit':   phi_i_doit,
            '_latex': phi_i_latex if not ij else phi_ij_latex
        }
    )


def is_lower_matrix(M: sp.Matrix):
    return sum([
        sp.Abs(M[i, j]) for i, j in zip(*np.triu_indices(M.cols, 1))
    ]) == 0


def is_strictly_lower_matrix(M: sp.Matrix):
    return sum([
        sp.Abs(M[i, j]) for i, j in zip(*np.triu_indices(M.cols))
    ]) == 0


def label_to_id(label: str):
    r = label.lower()
    replacements = [
        (" ", "_"),
        ("(", ""),
        (")", ""),
        (",", ""),
        ("-", ""),
        ("/", ""),
    ]
    for old, new in replacements:
        r = r.replace(old, new)
    return r


def tag(butcher: dict):
    if is_strictly_lower_matrix(butcher['A']):
        return "eRK"
    elif is_lower_matrix(butcher['A']):
        return "diRK"
    else:
        return "iRK"


def extract_method(file_list: str):
    for filename in file_list:
        with open(filename, 'r') as f:
            data = json.load(f)
        butcher = {
            'label': data['label'],
            'A': sp.Matrix([
                [sp.parse_expr(str(aij), transformations="all") for aij in ai]
                for ai in data['A']
            ]),
            'b': sp.Matrix([
                sp.parse_expr(str(bi), transformations="all") for bi in data['b']
            ]),
            'c': sp.Matrix([
                sp.parse_expr(str(ci), transformations="all") for ci in data['c']
            ])
        }
        if 'b2' in data:
            butcher['b2'] = sp.Matrix([
                sp.parse_expr(str(bi), transformations="all") for bi in data['b2']
            ])
        butcher['tag'] = data['tag'] if 'tag' in data else tag(butcher)

        yield butcher


class multisplit_list:
    """
      functor to split into multiple sublist a list

      condition should return a int between 0 and n-1
    """

    def __init__(self, n: int):
        self.n = n

    def __call__(self, seq: list, condition):
        r = [[] for _ in range(self.n)]
        for x in seq:
            r[condition(x)].append(x)
        return r


def tag_id(rk):
    select = {
        'eRK': 0,
        'expRK': 1,
        'diRK': 2,
        'iRK': 3
    }
    return select[rk['tag']]


def expRK_code_skeleton(X: list, c: list):
    z = sp.symbols('z')

    r = {
        'type': [],
        'code': []
    }

    for x in X:
        for symbol in x.free_symbols:
            # substitute all phi function
            x = x.subs(symbol, phi(
                *map(int, str(symbol).split("_")[1:]), c=c)(z))

        if len(x.free_symbols) == 0:
            r['type'].append("value_t")
            r['code'].append(x.evalf())
        else:
            r['type'].append("func_t")
            r['code'].append("[](linear_t && {}) -> linear_t {{ return {}; }}".format(
                str(z),
                sp.cxxcode(x.doit().simplify())
            ))

    return r


def prepare_expRK(rk: dict, Ndigit: int):
    print(rk['label']+" "*10, end="\r")
    extra_tableaus = ['b']
    if 'b2' in rk:
        extra_tableaus.append('b2')

    r = {
        x: expRK_code_skeleton(rk[x], rk['c'])
        for x in extra_tableaus
    }
    r['A'] = expRK_code_skeleton([rk['A'][i, j] for i, j in zip(
        *np.tril_indices(rk['A'].cols, k=-1))], rk['c'])

    c = rk['c'].evalf(n=Ndigit).T.tolist()[0]
    r['c'] = " , ".join(map(str, c))

    r['label'] = rk['label']
    r['id'] = label_to_id(r['label'])

    return r


def prepare_RK(rk: dict, Ndigit: int):
    print(rk['label']+" "*10, end="\r")
    sys.path.append("../analysis")
    from analysis import rk_butcher

    tmp_rk = rk_butcher(
        label=rk['label'],
        A=rk['A'], b=rk['b'], c=rk['c']
    )

    A = rk['A'].evalf(n=Ndigit).tolist()
    b = rk['b'].evalf(n=Ndigit).T.tolist()[0]
    c = rk['c'].evalf(n=Ndigit).T.tolist()[0]

    r = {
        'label': rk['label'],
        'id': tmp_rk.id,
        'A': [" , ".join(map(str, ai)) for ai in A],
        'b': " , ".join(map(str, b)),
        'c': " , ".join(map(str, c)),
        'order': tmp_rk.order
    }

    if 'b2' in rk:
        b2 = rk['b2'].evalf(n=Ndigit).T.tolist()[0]
        r['b2'] = " , ".join(map(str, b2))

    return r


def sformat(value, fmt, attribute=None):
    """
      filter for Jinja2 to transform a list into a list of string with format `fmt`
    """
    if attribute is not None:
        for elm in value:
            yield (fmt.format(elm[attribute]))
    else:
        for elm in value:
            yield (fmt.format(elm))


parser = argparse.ArgumentParser(
    description="code generator of Runge-Kutta method from them Butcher tableau")
parser.add_argument('FILE', nargs='*',
                    help="Files which contain a Butcher tableau")
parser.add_argument('-o', '--output', type=str,
                    help="Name of output file header C++")
parser.add_argument('--Ndigit', type=int, default=36, required=False,
                    help="number of digit in output Butcher tableau [default: 36]")
parser.add_argument('-d', '--doc', required=False, action='store_true')
parser.add_argument('-do', '--doc-output', type=str, required=False,
                    help="Name of output RST file")

if __name__ == '__main__':
    args = parser.parse_args()

    list_erk, list_exprk, list_dirk, list_irk = multisplit_list(
        4)(extract_method(args.FILE), tag_id)

    cg_list_erk = [prepare_RK(rk, args.Ndigit) for rk in list_erk]
    cg_list_exprk = [prepare_expRK(rk, args.Ndigit) for rk in list_exprk]
    cg_list_dirk = [prepare_RK(rk, args.Ndigit) for rk in list_dirk]

    print("code generation")

    local_dir = os.path.dirname(os.path.abspath(__file__))

    import jinja2
    jinja2.filters.FILTERS['sformat'] = sformat
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(local_dir))
    template = env.get_template("template/tpl_butcher_methods.hxx")

    with open(args.output, 'w') as butcher_hxx:
        butcher_hxx.write(
            template.render(list_erk=cg_list_erk,
                            list_exprk=cg_list_exprk, list_dirk=cg_list_dirk)
        )

    if args.doc:
        template_doc = env.get_template("template/tpl_doc_algorithms.rst")

        with open(args.doc_output, 'w') as file:
            file.write(
                template_doc.render(
                    list_erk=cg_list_erk, list_exprk=cg_list_exprk, list_dirk=cg_list_dirk)
            )
