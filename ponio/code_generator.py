# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#! /usr/bin/env python

import json
import argparse
import os

import sympy as sp
import numpy as np


class butcher_tableau:
    """
        class to represent a Butcher tableau and miscellaneous information about the method

        storing values:
        - A: matrix of Butcher tableau
        - b: vector of the last line of Butcher tableau
        - b2: optional last line of Butcher tableau for embedded method
        - c: vector of time coefficients of the Butcher tableau
        - tag: information about the type of method (for code generation)
        - doi: optional information about bibliography
    """

    def __init__(self, label, A, b, c, b2=None, tag=None, doi=None, **kwargs):
        self.label = label
        self.A = sp.Matrix([
            self._parse_vector(ai)
            for ai in A
        ])
        self.b = sp.Matrix(self._parse_vector(b))
        self.c = sp.Matrix(self._parse_vector(c))
        self.is_embedded = b2 is not None
        self.b2 = sp.Matrix(self._parse_vector(
            b2)) if self.is_embedded else None
        self.tag = tag
        self.doi = doi
        self.vphantom = ""

    @classmethod
    def from_file(cls, filename: str):
        return cls(**json.load(open(filename)))

    @classmethod
    def from_json(cls, json_dict):
        return cls(**json_dict)

    @classmethod
    def _parse_vector(cls, x):
        """
        Helper function to parse a vector of string into SymPy expression
        """
        return [sp.parse_expr(str(xi)) for xi in x]

    @property
    def id(self):
        """
        Returns an id from label
        """

        r = self.label.lower()
        replacements = [
            (" ", "_"),
            ("(", ""),
            (")", ""),
            ("[", ""),
            ("]", ""),
            (",", ""),
            ("-", ""),
            ("/", ""),
        ]
        for old, new in replacements:
            r = r.replace(old, new)
        return r

    @classmethod
    def line_repr_latex(cls, line, func=sp.latex, vphantom=""):
        """
        Helper to display a line in LaTeX

        :param line: vector to display in line
        :param func: function to call on each element of `line`
        :param vphantom: optional parameter to add a `vphantom` in LaTeX representation
        """
        return vphantom + " & ".join(map(func, line))

    def _repr_latex_(self, **kwargs):
        """
        Representation in notebook
        """

        format = "{{ c | {} }}".format('c'*len(self.A.row(0)))
        cA = [
            "{ci} & {aij} \\\\ \n".format(
                ci=sp.latex(ci),
                aij=self.line_repr_latex(self.A.row(i), lambda a_ij: sp.latex(
                    a_ij) if a_ij != 0 else " ", vphantom=self.vphantom)
            )
            for i, ci in enumerate(self.c)
        ]

        bb2 = self.line_repr_latex(self.b, vphantom=self.vphantom)

        if self.is_embedded:
            bb2 += f"\\\\ \n \\hline & {self.line_repr_latex(self.b2, vphantom=self.vphantom)}"

        return f"""\\begin{{array}}{format}
            {' '.join(cA)}\\hline
                & {bb2}
        \\end{{array}}"""

    def __iter__(self):
        yield 'label', self.label
        yield 'id', self.id

        if self.tag == "expRK":
            yield 'A', expRK_code_skeleton([self.A[i, j] for i, j in zip(
                *np.tril_indices(self.A.cols, k=-1))], self.c)
            yield 'b', expRK_code_skeleton(self.b, self.c)
            if self.is_embedded:
                yield 'b2', expRK_code_skeleton(self.b2, self.c)
        else:
            yield 'A', self.A.evalf().tolist()
            yield 'b', self.b.T.evalf().tolist()[0]
            if self.is_embedded:
                yield 'b2', self.b2.T.evalf().tolist()[0]

        yield 'c', self.c.T.evalf().tolist()[0]

        yield 'is_embedded', self.is_embedded

        if self.tag == "expRK":
            yield 'butcher', {
                'A': [[sp.latex(aij).replace("phi", "varphi") for aij in ai] for ai in self.A.tolist()],
                'b': [sp.latex(bi).replace("phi", "varphi") for bi in self.b.T.tolist()[0]],
                'c': list(map(sp.latex, self.c.T.tolist()[0]))
            }
        else:
            yield 'butcher', {
                'A': [
                    list(map(sp.latex, ai)) for ai in self.A.tolist()
                ],
                'c': list(map(sp.latex, self.c.T.tolist()[0]))
            } | {
                b: list(map(sp.latex, getattr(self, b).T.tolist()[0])) for b in ["b", "b2"] if getattr(self, b) is not None
            }

        if self.doi is not None:
            yield 'bib', doi_bib(self.doi)

        computer_order = get_computer_order(self.tag)
        R = computer_order.stability_function(self)
        yield 'stability_function', sp.latex(R(*R.signature))

        if not hasattr(self, 'order'):
            self.order = computer_order.order(self)
        yield 'order', self.order


class pair_butcher_tableau:
    """
    class to represent a pair of Butcher tableaus for additive Runge-Kutta methods
    """

    def __init__(self, label, *, explicit, implicit, tag=None, doi=None, **kwargs):
        self.label = label
        self.explicit = butcher_tableau(
            label=f"{label}-ex", **explicit, tag="eRK")
        self.implicit = butcher_tableau(
            label=f"{label}-im", **implicit, tag="diRK")
        self.is_embedded = self.explicit.b2 is not None and self.implicit.b2 is not None
        self.tag = tag
        self.doi = doi

        self.explicit.vphantom = r"\vphantom{ \frac{\sqrt{1}}{\sqrt{1}} }"
        self.implicit.vphantom = r"\vphantom{ \frac{\sqrt{1}}{\sqrt{1}} }"

    @classmethod
    def from_json(cls, json_dict):
        return cls(**json_dict)

    def _repr_latex_(self, **kwargs):
        return f"{self.explicit._repr_latex_()}\\quad{self.implicit._repr_latex_()}"

    @property
    def id(self):
        r = self.label.lower()
        replacements = [
            (" ", "_"),
            ("(", ""),
            (")", ""),
            ("[", ""),
            ("]", ""),
            (",", ""),
            ("-", ""),
            ("/", ""),
        ]
        for old, new in replacements:
            r = r.replace(old, new)
        return r

    def __iter__(self):
        yield 'label', self.label
        yield 'id', self.id
        yield 'explicit', dict(self.explicit)
        yield 'implicit', dict(self.implicit)

        yield 'is_embedded', self.is_embedded

        if self.doi is not None:
            yield 'bib', doi_bib(self.doi)

        computer_order = get_computer_order(self.tag)
        R = computer_order.stability_function(self)
        yield 'stability_function', sp.latex(R(*R.signature))

        if not hasattr(self, 'order'):
            self.order = computer_order.order(self)
        yield 'order', self.order


tags = ['eRK', 'expRK', 'diRK', 'iRK', 'aRK']


class rk_order:
    @classmethod
    def stability_function_Ab(cls, A, b):
        I = sp.eye(*A.shape)
        ones = sp.ones(*b.shape)

        z = sp.Dummy("z")
        R = 1 + z*(b.T*(I-z*A).inv()*ones)[0, 0]

        return sp.Lambda(z, R)

    @classmethod
    def stability_function(cls, rk):
        return cls.stability_function_Ab(rk.A, rk.b)

    @classmethod
    def order_star_function_Ab(cls, A, b):
        R = cls.stability_function_Ab(A, b)
        z, = R.signature

        return sp.Lambda(z, sp.Abs(R(z)*sp.exp(-z)))

    @classmethod
    def order_star_function(cls, rk):
        return cls.order_star_function_Ab(rk.A, rk.b)

    @classmethod
    def order_Ab(cls, A, b):
        order_star_func = cls.order_star_function_Ab(A, b)
        z, = order_star_func.signature

        theta = np.linspace(0., 2.*np.pi, 50)
        rho = 0.1

        order_star_eval = sp.lambdify(
            z, order_star_func(z)-1, 'numpy')(rho*np.exp(1j*theta))
        return sum(order_star_eval[1:]*order_star_eval[:-1] < 0.)//2 - 1

    @classmethod
    def order(cls, rk):
        return cls.order_Ab(rk.A, rk.b)


class ark_order:
    @classmethod
    def stability_function_Ab(cls, A_ex, b_ex, A_im, b_im):
        I = sp.eye(*A_ex.shape)
        ones = sp.ones(*b_ex.shape)

        z_ex, z_im = sp.Dummy("z_e"), sp.Dummy("z_i")
        R = 1 + ((z_ex*b_ex.T + z_im*b_im.T) *
                 (I - z_ex*A_ex - z_im*A_im).inv()*ones)[0, 0]

        return sp.Lambda((z_ex, z_im), R.expand())

    @classmethod
    def stability_function(cls, ark):
        return cls.stability_function_Ab(ark.explicit.A, ark.explicit.b, ark.implicit.A, ark.implicit.b)

    @classmethod
    def order_Ab(cls, A_ex, b_ex, A_im, b_im):
        R = cls.stability_function_Ab(A_ex, b_ex, A_im, b_im)
        z_ex, z_im = R.signature

        max_order = max(rk_order.order_Ab(A_ex, b_ex),
                        rk_order.order_Ab(A_im, b_im)) + 2

        coeff = sp.Poly(
            sp.series(sp.exp(z_ex+z_im), z_ex+z_im, 0, max_order).removeO()
            -
            sp.series(R(z_ex, z_im), z_im, 0, max_order).removeO(),
            z_ex, z_im
        ).as_dict()
        # for a polynomial as sum( a_{j,k}ze^j zi^k, j,k ) coeff is (j,k)=>a_{j,k}

        # sort coeff by order of 2 variables polynomial
        coeff_order = sorted(coeff.items(), key=lambda item: sum(item[0]))

        # get the first item where coeff is not zero (numerical zero)
        return next(sum(item[0]) for item in coeff_order if item[1] > 1e-12) - 1

    @classmethod
    def order(cls, ark):
        return cls.order_Ab(ark.explicit.A, ark.explicit.b, ark.implicit.A, ark.implicit.b)


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
            return ((phi(i-1)(z) - phi(i-1)(0))/z).expand().simplify()

    def phi_i_latex(self, printer, **kwargs):
        return f"\\varphi_{{{i}}}"
        # return str( type(self).mro()[0] ).replace(" ","\ ")

    def phi_ij_latex(self, printer, **kwargs):
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


class exprk_order:
    @classmethod
    def stability_function_Ab(cls, A, b, c):
        z_l, z_n = sp.Dummy("z_l"), sp.Dummy("z_n")

        I = sp.eye(*A.shape)
        ones = sp.ones(*b.shape)

        _A = A.subs({
            symbol: phi(*map(int, str(symbol).split("_")[1:]), c=c)(z_l)
            for symbol in A.free_symbols
        })
        _b = b.subs({
            symbol: phi(*map(int, str(symbol).split("_")[1:]))(z_l)
            for symbol in b.free_symbols
        })

        R = ((1 + z_n*(_b.T*(I - z_n*_A).inv()*(I + z_l*_A)*ones)
             [0, 0] + z_l*(_b.T*ones)[0, 0]).expand()).collect(z_l).collect(z_n)

        return sp.Lambda((z_l, z_n), R)

    @classmethod
    def stability_function(cls, rk):
        return cls.stability_function_Ab(rk.A, rk.b, rk.c)

    @classmethod
    def order_Ab(cls, A, b, c):
        R = cls.stability_function_Ab(A, b, c)
        z_l, z_n = R.signature

        R_zl0 = sp.limit(R(z_l, z_n).doit().expand(), z_l, 0)
        max_order = sp.Poly(R_zl0, z_n).degree() + 2

        coeff = sp.Poly(
            sp.series(sp.exp(z_n), z_n, 0, max_order).removeO()
            -
            R_zl0
        ).as_dict()

        # sort coeff by order of 2 variables polynomial
        coeff_order = sorted(coeff.items(), key=lambda item: sum(item[0]))

        # get the first item where coeff is not zero (numerical zero)
        return next(sum(item[0]) for item in coeff_order if item[1] > 1e-12) - 1

    @classmethod
    def order(cls, rk):
        return cls.order_Ab(rk.A, rk.b, rk.c)


def get_computer_order(tag):
    if tag in ('eRK', 'diRK', 'iRK'):
        return rk_order
    if tag in ('aRK'):
        return ark_order
    if tag in ('expRK'):
        return exprk_order


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
            r['code'].append(
                f"[](linear_t && {z}) -> linear_t {{ return {sp.cxxcode(x.doit().simplify())}; }}")

    return r


def doi_bib_crossref(doi: str):
    """
        return bibliography reference from doi with a request to `api.crossref.org`
    """
    import urllib.request

    try:
        with urllib.request.urlopen(f"https://api.crossref.org/works/{doi}") as response:
            message = json.loads(response.read())['message']

            author = " & ".join(
                [f"{auth['family']}, {auth['given']}" for auth in message['author']])
            title = message['title'][0]
            pubdate = message['published']['date-parts'][0][0]
            publisher = message['short-container-title'][0] if len(
                message['short-container-title']) > 0 else message['publisher']

        return {
            'url': message['URL'],
            'bib': f"{author}, \"{title}\", in: {publisher} ({pubdate})"
        }
    except urllib.error.HTTPError as e:
        print(f"HTTPError: {e.code}")

        return {
            'url': f"https://www.doi.org/{doi}",
            'bib': f"{doi} [preprint]"
        }


def doi_bib_offline(doi: str):
    """
        return bibliography reference as doi and url to doi website
    """
    return {
        'url': f'https://doi.org/{doi}',
        'bib': doi
    }


# make a request to get bibliography by default
doi_bib = doi_bib_crossref


class multisplit_list:
    def __init__(self, tags, key):
        self.tags = tags
        self.key = key

    def __call__(self, seq):
        r = {tag: [] for tag in self.tags}
        for x in seq:
            r[self.key(x)].append(x)
        return r


def extract_method(file_list):
    for filename in file_list:
        with open(filename, 'r') as f:
            data = json.load(f)

        if data['tag'] in ('eRK', 'expRK', 'diRK', 'iRK'):
            rk = butcher_tableau.from_json(data)
        elif data['tag'] in ('aRK'):
            rk = pair_butcher_tableau.from_json(data)

        yield rk


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
parser.add_argument('-d', '--doc', required=False, action='store_true',
                    help="generate also Sphinx documentation (api/algorithm.rst file)")
parser.add_argument('-do', '--doc-output', type=str, required=False,
                    help="Name of output RST directory")
parser.add_argument('--offline', required=False, action='store_true',
                    help="make an offline generation (without bibliography research in Doxygen comments)")

if __name__ == '__main__':
    import jinja2

    args = parser.parse_args()

    if args.offline:
        doi_bib = doi_bib_offline

    all_meths = multisplit_list(tags, lambda rk: rk.tag)(
        extract_method(args.FILE))
    all_meths.pop('iRK', None)  # not implemented method in ponio

    total_meth = 0
    for tag, meths in all_meths.items():
        print(f"{tag:>5} {len(meths)}")
        total_meth += len(meths)
    print('  lRK', len(all_meths['eRK']), "(from eRK methods)")
    total_meth += len(all_meths['eRK'])
    print("-"*8)
    print("total", total_meth)

    jinja2.filters.FILTERS['sformat'] = sformat
    local_dir = os.path.dirname(os.path.abspath(__file__))
    template_dir = os.path.join(local_dir, "template")
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(template_dir)
    )
    template = env.get_template("tpl_butcher_methods.cpp.jinja2")

    def dict_and_log(rk):
        print(f"{rk.label:>30}", end="\r")
        return dict(rk)

    list_erk = [dict_and_log(rk) for rk in all_meths['eRK']]
    list_dirk = [dict_and_log(rk) for rk in all_meths['diRK']]
    list_exprk = [dict_and_log(rk) for rk in all_meths['expRK']]
    list_ark = [dict_and_log(rk) for rk in all_meths['aRK']]

    with open(args.output, 'w') as butcher_hxx:
        butcher_hxx.write(template.render(
            list_erk=list_erk,
            list_dirk=list_dirk,
            list_exprk=list_exprk,
            list_ark=list_ark
        ))

    if args.doc:
        # main list
        template_doc = env.get_template("tpl_doc_algorithms.rst")

        with open(f"{args.doc_output}/list_algorithm.rst", 'w') as file:
            file.write(
                template_doc.render(
                    list_erk=list_erk,
                    list_dirk=list_dirk,
                    list_exprk=list_exprk,
                    list_ark=list_ark
                )
            )

        # sublists
        for tpl in ("erk", "dirk", "lrk", "dp", "exprk", "ark"):
            template_doc = env.get_template(f"tpl_doc_{tpl}.rst")

            with open(f"{args.doc_output}/list_alg_{tpl}.rst", 'w') as file:
                file.write(
                    template_doc.render(
                        list_erk=list_erk,
                        list_dirk=list_dirk,
                        list_exprk=list_exprk,
                        list_ark=list_ark
                    )
                )
