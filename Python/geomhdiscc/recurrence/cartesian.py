"""Module to compute the recurrence relations for the spectral operators of the a cartesian direction"""

from __future__ import division
from __future__ import unicode_literals

import sympy
import geomhdiscc.recurrence.symbolic as symbolic


def i1():
    """Cartesian first integral operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Cartesian second integral operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2lapl():
    """Cartesian second integral of laplacian operator"""

    # Setup terms in recurrence
    k, l = sympy.symbols('k l')
    terms = [{'q':2, 'p':0, 'd':2, 'c':1},{'q':2, 'p':0, 'd':0, 'c':-k**2},{'q':2, 'p':0, 'd':0, 'c':-l**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2laplh():
    """Cartesian second integral of horizontal laplacian operator"""

    # Setup terms in recurrence
    k = sympy.Symbol('k')
    terms = [{'q':2, 'p':0, 'd':2, 'c':1},{'q':2, 'p':0, 'd':0, 'c':-k**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4():
    """Cartesian fourth integral operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4d2():
    """Cartesian fourth integral of second derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':2, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl():
    """Cartesian fourth integral of laplacian operator"""

    # Setup terms in recurrence
    k, l = sympy.symbols('k l')
    terms = [{'q':4, 'p':0, 'd':2, 'c':1},{'q':4, 'p':0, 'd':0, 'c':-k**2},{'q':4, 'p':0, 'd':0, 'c':-l**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4laplh():
    """Cartesian fourth integral of horizontal laplacian operator"""

    # Setup terms in recurrence
    k = sympy.Symbol('k')
    terms = [{'q':4, 'p':0, 'd':2, 'c':1},{'q':4, 'p':0, 'd':0, 'c':-k**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl2():
    """Cartesian fourth integral of bilaplacian operator"""

    # Setup terms in recurrence
    k, l = sympy.symbols('k l')
    terms = [{'q':4, 'p':0, 'd':4, 'c':1},{'q':4, 'p':0, 'd':0, 'c':k**4},{'q':4, 'p':0, 'd':0, 'c':l**4},{'q':4, 'p':0, 'd':2, 'c':-2*k**2},{'q':4, 'p':0, 'd':2, 'c':-2*l**2},{'q':4, 'p':0, 'd':0, 'c':2*k**2*l**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl2h():
    """Cartesian fourth integral of horizontal bilaplacian operator"""

    # Setup terms in recurrence
    k = sympy.Symbol('k')
    terms = [{'q':4, 'p':0, 'd':4, 'c':1},{'q':4, 'p':0, 'd':0, 'c':k**4},{'q':4, 'p':0, 'd':2, 'c':-2*k**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
