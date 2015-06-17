"""Module to compute the recurrence relations for the spectral Jacobi operators of the radial direction in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import geomhdiscc.recurrence.symbolic_jacobi as symbolic

#a = sympy.Symbol('a')
#b = sympy.Symbol('b')
l = sympy.Symbol('l')
n = sympy.Symbol('n')

def x1():
    """Sphere x operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r2():
    """Sphere r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':sympy.Rational(1,2)},{'q':0, 'p':0, 'd':0, 'c':sympy.Rational(1,2)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r4():
    """Sphere r^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':2, 'd':0, 'c':sympy.Rational(1,4)},{'q':0, 'p':1, 'd':0, 'c':sympy.Rational(1,2)},{'q':0, 'p':0, 'd':0, 'c':sympy.Rational(1,4)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1():
    """Sphere i1 operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':sympy.Rational(1,4)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Sphere i2 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':sympy.Rational(1,16)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4():
    """Sphere i4 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':sympy.Rational(1,4**4)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2lapl():
    """Sphere i2lapl operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':1, 'd':2, 'c':sympy.Rational(1,2)}, {'q':2, 'p':0, 'd':2, 'c':sympy.Rational(1,2)}, {'q':2, 'p':0, 'd':1, 'c':l*sympy.Rational(1,2) + sympy.Rational(3,4)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl():
    """Sphere i2lapl operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':1, 'd':2, 'c':sympy.Rational(1,32)}, {'q':4, 'p':0, 'd':2, 'c':sympy.Rational(1,32)}, {'q':4, 'p':0, 'd':1, 'c':l*sympy.Rational(1,32) + sympy.Rational(3,64)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl2():
    """Sphere i2lapl2 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':2, 'd':4, 'c':sympy.Rational(1,4)}, {'q':4, 'p':1, 'd':4, 'c':sympy.Rational(1,2)}, {'q':4, 'p':0, 'd':4, 'c':sympy.Rational(1,4)}, {'q':4, 'p':1, 'd':3, 'c':l*sympy.Rational(1,2) + sympy.Rational(5,4)}, {'q':4, 'p':0, 'd':3, 'c':l*sympy.Rational(1,2) + sympy.Rational(5,4)}, {'q':4, 'p':0, 'd':2, 'c':l**2*sympy.Rational(1,4) + l + sympy.Rational(15,16)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
