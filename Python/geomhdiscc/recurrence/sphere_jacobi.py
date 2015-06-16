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
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Sphere i2 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4():
    """Sphere i4 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2lapl():
    """Sphere i2lapl operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':1, 'd':2, 'c':8}, {'q':2, 'p':0, 'd':2, 'c':8}, {'q':2, 'p':0, 'd':1, 'c':4*(2*l+3)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl():
    """Sphere i2lapl operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':1, 'd':2, 'c':8}, {'q':4, 'p':0, 'd':2, 'c':8}, {'q':4, 'p':0, 'd':1, 'c':4*(2*l+3)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl2():
    """Sphere i2lapl2 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':2, 'd':4, 'c':64}, {'q':4, 'p':1, 'd':4, 'c':128}, {'q':4, 'p':0, 'd':4, 'c':64}, {'q':4, 'p':1, 'd':3, 'c':64*(2*l+5)}, {'q':4, 'p':0, 'd':3, 'c':64*(2*l+5)}, {'q':4, 'p':0, 'd':2, 'c':16*(4*l**2 + 16*l + 15)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
