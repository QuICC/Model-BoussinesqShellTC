"""Module to compute the recurrence relations for the spectral Jacobi operators of the radial direction in a cylinder"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import geomhdiscc.recurrence.symbolic_jacobi as mod

m = sympy.Symbol('m')
n = sympy.Symbol('n')

symbolic = mod.SymbolicJacobi(a = -sympy.Rational(1,2), b = m - sympy.Rational(1,2))

def x1():
    """Cylinder x operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r2():
    """Cylinder r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':sympy.Rational(1,2)},{'q':0, 'p':0, 'd':0, 'c':sympy.Rational(1,2)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r4():
    """Cylinder r^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':2, 'd':0, 'c':sympy.Rational(1,4)},{'q':0, 'p':1, 'd':0, 'c':sympy.Rational(1,2)},{'q':0, 'p':0, 'd':0, 'c':sympy.Rational(1,4)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1():
    """Cylinder i1 (i1r1) operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Cylinder i2 ((i1r1)^2) operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2laplh():
    """Cylinder i2laplh ((i1r1)^2 laplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':2, 'p':1, 'd':2, 'c':8},
            {'q':2, 'p':0, 'd':2, 'c':8},
            {'q':2, 'p':0, 'd':1, 'c':8*(m + 1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4():
    """Cylinder i4 ((i1r1)^4) operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4laplh():
    """Cylinder i4laplh ((i1r1)^4 laplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':4, 'p':1, 'd':2, 'c':8},
            {'q':4, 'p':0, 'd':2, 'c':8},
            {'q':4, 'p':0, 'd':1, 'c':8*(m + 1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl2h():
    """Cylinder i4lapl2h ((i1r1)^4 bilaplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':4, 'p':2, 'd':4, 'c':64},
            {'q':4, 'p':1, 'd':4, 'c':128},
            {'q':4, 'p':0, 'd':4, 'c':64},
            {'q':4, 'p':1, 'd':3, 'c':128*(m+2)},
            {'q':4, 'p':0, 'd':3, 'c':128*(m+2)},
            {'q':4, 'p':0, 'd':2, 'c':64*(m+2)*(m+1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6():
    """Cylinder i6 ((i1r1)^6) operator"""

    # Setup terms in recurrence
    terms = [{'q':6, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6laplh():
    """Cylinder i6laplh ((i1r1)^4 laplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':6, 'p':1, 'd':2, 'c':8},
            {'q':6, 'p':0, 'd':2, 'c':8},
            {'q':6, 'p':0, 'd':1, 'c':8*(m + 1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6lapl2h():
    """Cylinder i6lapl2 ((i1r1)^6 bilaplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':6, 'p':2, 'd':4, 'c':64},
            {'q':6, 'p':1, 'd':4, 'c':128},
            {'q':6, 'p':0, 'd':4, 'c':64},
            {'q':6, 'p':1, 'd':3, 'c':128*(m+2)},
            {'q':6, 'p':0, 'd':3, 'c':128*(m+2)},
            {'q':6, 'p':0, 'd':2, 'c':64*(m+2)*(m+1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6lapl3h():
    """Cylinder i6lapl3 ((i1r1)^6 trilaplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':6, 'p':3, 'd':6, 'c':512},
            {'q':6, 'p':2, 'd':6, 'c':1536},
            {'q':6, 'p':1, 'd':6, 'c':1536},
            {'q':6, 'p':0, 'd':6, 'c':512},
            {'q':6, 'p':2, 'd':5, 'c':1536*(m+3)},
            {'q':6, 'p':1, 'd':5, 'c':3072*(m+3)},
            {'q':6, 'p':0, 'd':5, 'c':1536*(m+3)},
            {'q':6, 'p':1, 'd':4, 'c':1536*(m+3)*(m+2)},
            {'q':6, 'p':0, 'd':4, 'c':1536*(m+3)*(m+2)},
            {'q':6, 'p':0, 'd':3, 'c':512*(m+3)*(m+2)*(m+1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
