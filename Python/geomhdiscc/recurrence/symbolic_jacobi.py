"""Module used for symbolic manipulations of sparse Jacobi expansions"""

from __future__ import division
from __future__ import unicode_literals

import sympy
import copy

l = sympy.Symbol('l')
n = sympy.Symbol('n')
a = -sympy.Rational(1,2)
b = l - sympy.Rational(1,2)

import geomhdiscc.recurrence.symbolic_base as base

class SymbolicJacobi(base.SymbolicBase):

    def norm(self, i):
        if i == 0:
            val = 1
        elif i > 0:
            val = (2*n + a + b + 1)
            for j in range(0, i):
                val = val*(n + j + a + 1)*(n + j + b + 1)/((2*n + 2*j + a + b + 3)*(n + j + a + b + 1)*(n + j + 1))
            val = sympy.sqrt(val.simplify())
        else:
            val = 1

        return val

    def spectral_monomial(self, p, f, asrow = True):
        # Some assertion for safety
        assert (p >= 0)

        # End recurrence: monomial order zero is identity
        if p == 0:
            recurrence = dict(f)
        else:
            prev = self.spectral_monomial(p-1, f, False)
            recurrence = dict()

            def cA(i):
                """A coefficient: x J_n =  A_n J_{n-1} - B_n J_{n} + C_n J_{n+1}"""
                
                num = 2*(n + i + a)*(n + i + b)
                den = (2*n + 2*i + a + b + 1)*(2*n + 2*i + a + b)
                return num/den

            def cB(i):
                """B coefficient: x J_n =  A_n J_{n-1} - B_n J_{n} + C_n J_{n+1}"""
               
                num = (a**2 - b**2)
                den = (2*n + 2*i + a + b + 2)*(2*n + 2*i + a + b)
                return num/den

            def cC(i):
                """C coefficient: x J_n =  A_n J_{n-1} - B_n J_{n} + C_n J_{n+1}"""

                num = 2*(n + i + 1)*(n + i + a + b + 1)
                den = (2*n + 2*i + a + b + 1)*(2*n + 2*i + a + b + 2)
                return num/den

            for i in prev.keys():
                recurrence[i-1] = recurrence.get(i-1,0) + cA(i)*prev[i]
                recurrence[i] = recurrence.get(i,0) - cB(i)*prev[i]
                recurrence[i+1] = recurrence.get(i+1,0) + cC(i)*prev[i]

            # Convert to row recurrence relation
            if asrow:
                old = dict(recurrence)
                recurrence = dict()
                for i in old.keys():
                    recurrence[-i] = old[i].subs(n, n-i)

            for i in recurrence.keys():
                recurrence[i] = recurrence[i].simplify().factor()
        return recurrence

    def spectral_integral(self, q, f, asrow = True):
        # Some assertion for safety
        assert (q >= 0)

        # End recurrence: integration "power" zero is identity
        if q == 0:
            recurrence = dict(f)
        else:
            prev = self.spectral_integral(q-1, f, False)
            recurrence = dict()

            def cA(i):
                """A coefficient: \int J_n =  A_n J_{n-1} + B_n J_{n} + C_n J_{n+1}"""
                
                num = -2*(n + i + a)*(n + i + b)
                den = (n + i + a + b)*(2*n + 2*i + a + b)*(2*n + 2*i + a + b + 1)
                return num/den

            def cB(i):
                """B coefficient: \int J_n =  A_n J_{n-1} + B_n J_{n} + C_n J_{n+1}"""
               
                num = 2*(a - b)
                den = (2*n + 2*i + a + b + 2)*(2*n + 2*i + a + b)
                return num/den

            def cC(i):
                """C coefficient: \int J_n =  A_n J_{n-1} + B_n J_{n} + C_n J_{n+1}"""

                num = 2*(n + i + a + b + 1)
                den = (2*n + 2*i + a + b + 1)*(2*n + 2*i + a + b + 2)
                return num/den

            for i in prev.keys():
                recurrence[i-1] = recurrence.get(i-1,0) + cA(i)*prev[i]
                recurrence[i] = recurrence.get(i,0) + cB(i)*prev[i]
                recurrence[i+1] = recurrence.get(i+1,0) + cC(i)*prev[i]

            # Convert to row recurrence relation
            if asrow:
                old = dict(recurrence)
                recurrence = dict()
                for i in old.keys():
                    recurrence[-i] = old[i].subs(n, n-i)

            for i in recurrence.keys():
                recurrence[i] = recurrence[i].simplify().factor()
        return recurrence
