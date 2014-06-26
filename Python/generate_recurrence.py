from __future__ import division
from __future__ import unicode_literals

import sympy
import symbolic

print("################################################")
print("################ Cartesian 1D ##################")
print("################################################")

print("i1:")
terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i2:")
terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i2lapl:")
k, l = sympy.symbols('k l')
terms = [{'q':2, 'p':0, 'd':2, 'c':1},{'q':2, 'p':0, 'd':0, 'c':-k**2},{'q':2, 'p':0, 'd':0, 'c':-l**2}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i2laplh:")
k = sympy.Symbol('k')
terms = [{'q':2, 'p':0, 'd':2, 'c':1},{'q':2, 'p':0, 'd':0, 'c':-k**2}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4:")
terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4d2:")
terms = [{'q':4, 'p':0, 'd':2, 'c':1}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4lapl:")
k, l = sympy.symbols('k l')
terms = [{'q':4, 'p':0, 'd':2, 'c':1},{'q':4, 'p':0, 'd':0, 'c':-k**2},{'q':4, 'p':0, 'd':0, 'c':-l**2}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4laplh:")
k = sympy.Symbol('k')
terms = [{'q':4, 'p':0, 'd':2, 'c':1},{'q':4, 'p':0, 'd':0, 'c':-k**2}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4lapl2:")
k, l = sympy.symbols('k l')
terms = [{'q':4, 'p':0, 'd':4, 'c':1},{'q':4, 'p':0, 'd':0, 'c':k**4},{'q':4, 'p':0, 'd':0, 'c':l**4},{'q':4, 'p':0, 'd':2, 'c':-2*k**2},{'q':4, 'p':0, 'd':2, 'c':-2*l**2},{'q':4, 'p':0, 'd':0, 'c':2*k**2*l**2}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4lapl2h:")
k = sympy.Symbol('k')
terms = [{'q':4, 'p':0, 'd':4, 'c':1},{'q':4, 'p':0, 'd':0, 'c':k**4},{'q':4, 'p':0, 'd':2, 'c':-2*k**2}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("################################################")
print("################ Spherical shell ###############")
print("################################################")
print("i2x2:")
terms = [{'q':2, 'p':2, 'd':0, 'c':1}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i2x2lapl:")
l = sympy.Symbol('l')
terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':2}, {'q':2, 'p':0, 'd':0, 'c':-l*(l+1)}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4x4:")
terms = [{'q':4, 'p':4, 'd':0, 'c':1}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4x4lapl:")
l = sympy.Symbol('l')
terms = [{'q':4, 'p':4, 'd':2, 'c':1}, {'q':4, 'p':3, 'd':1, 'c':2}, {'q':4, 'p':2, 'd':0, 'c':-l*(l+1)}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4x4lapl2:")
l = sympy.Symbol('l')
terms = [{'q':4, 'p':4, 'd':4, 'c':1}, {'q':4, 'p':3, 'd':3, 'c':4}, {'q':4, 'p':2, 'd':2, 'c':-2*l*(l+1)}, {'q':4, 'p':0, 'd':0, 'c':(l-1)*l*(l+1)*(l+2)}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("################################################")
print("#################### Sphere ####################")
print("################################################")
print("i2x2:")
terms = [{'q':2, 'p':2, 'd':0, 'c':1}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")

print("i2x2lapl:")
l = sympy.Symbol('l')
terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':2}, {'q':2, 'p':0, 'd':0, 'c':-l*(l+1)}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")

print("i4x4:")
terms = [{'q':4, 'p':4, 'd':0, 'c':1}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")

print("i4x4lapl:")
l = sympy.Symbol('l')
terms = [{'q':4, 'p':4, 'd':2, 'c':1}, {'q':4, 'p':3, 'd':1, 'c':2}, {'q':4, 'p':2, 'd':0, 'c':-l*(l+1)}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")

print("i4x4lapl2:")
l = sympy.Symbol('l')
terms = [{'q':4, 'p':4, 'd':4, 'c':1}, {'q':4, 'p':3, 'd':3, 'c':4}, {'q':4, 'p':2, 'd':2, 'c':-2*l*(l+1)}, {'q':4, 'p':0, 'd':0, 'c':(l-1)*l*(l+1)*(l+2)}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")

print("################################################")
print("################### Annulus ####################")
print("################################################")
print("i2x2:")
terms = [{'q':2, 'p':2, 'd':0, 'c':1}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i2x2lapl:")
m = sympy.Symbol('m')
terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':1}, {'q':2, 'p':0, 'd':0, 'c':-m*m}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4x4:")
terms = [{'q':4, 'p':4, 'd':0, 'c':1}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4x4lapl:")
m = sympy.Symbol('m')
terms = [{'q':4, 'p':4, 'd':2, 'c':1}, {'q':4, 'p':3, 'd':1, 'c':1}, {'q':4, 'p':2, 'd':0, 'c':-m**2}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("i4x4lapl2:")
m = sympy.Symbol('m')
terms = [{'q':4, 'p':4, 'd':4, 'c':1}, {'q':4, 'p':3, 'd':3, 'c':2}, {'q':4, 'p':2, 'd':2, 'c':-(1+2*m**2)}, {'q':4, 'p':1, 'd':1, 'c':(1+2*m**2)}, {'q':4, 'p':0, 'd':0, 'c':(m**2 - 4)*m**2}]
terms = symbolic.change_variable(terms, 'linear_r2x')
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")

print("################################################")
print("################### Cylinder  ##################")
print("################################################")
print("i2x2:")
terms = [{'q':2, 'p':2, 'd':0, 'c':1}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")

print("i2x2lapl:")
m = sympy.Symbol('m')
terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':1}, {'q':2, 'p':0, 'd':0, 'c':-m*m}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")

print("i4x4:")
terms = [{'q':4, 'p':4, 'd':0, 'c':1}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")

print("i4x4lapl:")
m = sympy.Symbol('m')
terms = [{'q':4, 'p':4, 'd':2, 'c':1}, {'q':4, 'p':3, 'd':1, 'c':1}, {'q':4, 'p':2, 'd':0, 'c':-m**2}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")

print("i4x4lapl2:")
m = sympy.Symbol('m')
terms = [{'q':4, 'p':4, 'd':4, 'c':1}, {'q':4, 'p':3, 'd':3, 'c':2}, {'q':4, 'p':2, 'd':2, 'c':-(1+2*m**2)}, {'q':4, 'p':1, 'd':1, 'c':(1+2*m**2)}, {'q':4, 'p':0, 'd':0, 'c':(m**2 - 4)*m**2}]
r = symbolic.build_recurrence(terms, {0:1})
n = sympy.Symbol('n')
print("   General:")
for k,rec in sorted(r.items()):
   print("\t" + str(k) + ": \t" + str(rec))
print("\n")
#print("   Even:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n)))
#print("\n")
#print("   Odd:")
#for k,rec in sorted(r.items()):
#   print("\t" + str(k/2) + ": \t" + str(rec.subs(n,2*n+1)))
#print("\n")
