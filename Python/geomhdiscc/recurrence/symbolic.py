"""Module used for symbolic manipulations of sparse chebyshev expansions"""

from __future__ import division
from __future__ import unicode_literals

import sympy
import copy

def integrate(term, base = None):
    # Some assertion for safety
    assert (term['q'] >= 0)
    assert (term['p'] >= 0)
    assert (term['d'] >= 0)
    assert (term['q'] >= term['d'])

    if base is None: base = dict()

    # End recurrence: derivative order is zero
    if term['d'] == 0:
        base[term['q'],term['p']] = base.get((term['q'],term['p']),0) + term['c']
    else:
        # Compute first part of integration by parts
        t = {'q':term['q'] - 1, 'p':term['p'], 'd':term['d'] - 1, 'c':term['c']}
        base = integrate(t, base)

        # Compute second part of integration by parts
        if term['p'] - 1 >= 0:
            t = {'q':term['q'], 'p':term['p'] - 1, 'd':term['d'] - 1, 'c':-term['p']*term['c']}
            base = integrate(t, base)
    return base


def spectral_monomial(p, f):
    # Some assertion for safety
    assert (p >= 0)

    # End recurrence: monomial order zero is identity
    if p == 0:
        recurrence = dict(f)
    else:
        prev = spectral_monomial(p-1, f)
        recurrence = dict()

        for i in prev.keys():
            recurrence[i-1] = recurrence.get(i-1,0) + sympy.Rational(1,2)*prev[i]
            recurrence[i+1] = recurrence.get(i+1,0) + sympy.Rational(1,2)*prev[i]
        for i in recurrence.keys():
            recurrence[i] = recurrence[i].simplify().factor()
    return recurrence


def spectral_integral(q, f, asrow = True):
    # Some assertion for safety
    assert (q >= 0)

    # End recurrence: integration "power" zero is identity
    if q == 0:
        recurrence = dict(f)
    else:
        prev = spectral_integral(q-1, f, False)
        recurrence = dict()

        n = sympy.Symbol('n')
        for i in prev.keys():
            recurrence[i-1] = recurrence.get(i-1,0) - (1/(2*(n+i-1)))*prev[i]
            recurrence[i+1] = recurrence.get(i+1,0) + (1/(2*(n+i+1)))*prev[i]

        # Convert to row recurrence relation
        if asrow:
            old = dict(recurrence)
            recurrence = dict()
            for i in old.keys():
                recurrence[-i] = old[i].subs(n, n-i)

        for i in recurrence.keys():
            recurrence[i] = recurrence[i].simplify().factor()
    return recurrence

def spectral_intgxmult(q, p, f):
    # Some assertion for safety
    assert (q >= 0)
    assert (p >= 0)

    # Compute monomial multiplication
    recurrence = spectral_monomial(p, f)

    # Compute integration
    recurrence = spectral_integral(q, recurrence)

    return recurrence


def change_variable(terms, type):
    new_terms = []
    if type == 'linear_r2x':
        a, b, x = sympy.symbols('a b x')
        for term in terms:
            poly = sympy.poly((a*x+b)**term['p'], x)
            pcoeffs = list(reversed(poly.coeffs()))
            for j,c in enumerate(pcoeffs):
                t = term.copy()
                t['c'] = c*t['c']*a**(t['q'] - t['d'])
                t['p'] = j
                new_terms.append(t)

    return new_terms


def build_recurrence(terms, fs):
    ops = dict();

    # Integrate and sum all terms
    for term in terms:
        op = integrate(term)
        for tup in op.keys():
            ops[tup] = ops.get(tup, 0) + op[tup]

    # Create recurrence relation
    recurrence = dict()
    for tup in ops.keys():
        for d,f in fs.items():
            rec = spectral_intgxmult(tup[0], tup[1], {d:f})
            for i in rec.keys():
                recurrence[i] = recurrence.get(i,0) + ops[tup]*rec[i]


    for i in recurrence.keys():
        recurrence[i] = recurrence[i].simplify().factor()
    return recurrence
