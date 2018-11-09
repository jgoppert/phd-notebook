from IPython.display import Latex
import sympy
import sympy.physics.mechanics as mech


def labelled_sympy_eq(expr, label):
    l_expr = mech.latex(expr)
    return Latex(r"""
\begin{{align}}
{l_expr:s}
\tag{{{label:s}}}\label{{eq:{label:s}}}
\end{{align}}""".format(**locals()))


def custom_latex_labelled_eq(lhs, expr):
    return Latex('$'+lhs+'='+mech.mlatex(expr)+'$')


def collect_within_terms(expr, f):
    new = 0
    for term in expr.as_ordered_terms():
        new += term.collect(f, sympy.factor)
    return new
