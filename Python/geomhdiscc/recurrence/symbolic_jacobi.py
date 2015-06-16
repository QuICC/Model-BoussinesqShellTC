"""Module used for symbolic manipulations of sparse Jacobi expansions"""

from __future__ import division
from __future__ import unicode_literals

import sympy
import copy

l = sympy.Symbol('l')
n = sympy.Symbol('n')
a = -sympy.Rational(1,2)
b = l - sympy.Rational(1,2)

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

def norm(i):
    if i == 0:
        norm = 1
    elif i > 0:
        norm = (2*n + a + b + 1)
        for j in range(0, i):
            norm = norm*(n + j + a + 1)*(n + j + b + 1)/((2*n + 2*j + a + b + 3)*(n + j + a + b + 1)*(n + j + 1))
        norm = sympy.sqrt(norm.simplify())
    else:
        norm = 1

    return norm

def spectral_monomial(p, f, asrow = True):
    # Some assertion for safety
    assert (p >= 0)

    # End recurrence: monomial order zero is identity
    if p == 0:
        recurrence = dict(f)
    else:
        prev = spectral_monomial(p-1, f, False)
        recurrence = dict()

        for i in prev.keys():
            recurrence[i-1] = recurrence.get(i-1,0) + 2*(n + i + a)*(n + i + b)/((2*n + 2*i + a + b + 1)*(2*n + 2*i + a + b))*prev[i]
            recurrence[i] = recurrence.get(i,0) - (a**2 - b**2)/((2*n + 2*i + a + b + 2)*(2*n + 2*i + a + b))*prev[i]
            recurrence[i+1] = recurrence.get(i+1,0) + 2*(n + i + 1)*(n + i + a + b + 1)/((2*n + 2*i + a + b + 1)*(2*n + 2*i + a + b + 2))*prev[i]

        # Convert to row recurrence relation
        if asrow:
            old = dict(recurrence)
            recurrence = dict()
            for i in old.keys():
                recurrence[-i] = old[i].subs(n, n-i)

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

        for i in prev.keys():
            recurrence[i-1] = recurrence.get(i-1,0) - 2*(n + i + a)*(n + i + b)/((n + i + a + b)*(2*n + 2*i + a + b)*(2*n + 2*i + a + b + 1))*prev[i]
            recurrence[i] = recurrence.get(i,0) + 2*(a - b)/((2*n + 2*i + a + b)*(2*n + 2*i + a + b + 2))*prev[i]
            recurrence[i+1] = recurrence.get(i+1,0) + 2*(n + i + a + b + 1)/((2*n + 2*i + a + b + 1)*(2*n + 2*i + a + b + 2))*prev[i]

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
