"""Module to compute the recurrence relations for the spectral operators of the radial direction in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy
import geomhdiscc.recurrence.symbolic as symbolic


def i2x2():
    """Sphere second integral of x^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':2, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
       print("\t" + str(k) + ": \t" + str(rec))
    print("\n")


def i2x2lapl():
    """Sphere second integral of x^2 laplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':2}, {'q':2, 'p':0, 'd':0, 'c':-l*(l+1)}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
       print("\t" + str(k) + ": \t" + str(rec))
    print("\n")


def i4x4():
    """Sphere fourth integral of x^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':4, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
       print("\t" + str(k) + ": \t" + str(rec))
    print("\n")


def i4x4lapl():
    """Sphere fourth integral of x^4 laplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':2, 'c':1}, {'q':4, 'p':3, 'd':1, 'c':2}, {'q':4, 'p':2, 'd':0, 'c':-l*(l+1)}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
       print("\t" + str(k) + ": \t" + str(rec))
    print("\n")


def i4x4lapl2():
    """Sphere fourth integral of x^4 bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':4, 'c':1}, {'q':4, 'p':3, 'd':3, 'c':4}, {'q':4, 'p':2, 'd':2, 'c':-2*l*(l+1)}, {'q':4, 'p':0, 'd':0, 'c':(l-1)*l*(l+1)*(l+2)}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
       print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
