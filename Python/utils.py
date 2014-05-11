"""Module provides generic functions for the sparse chebyshev representation."""

def build_diagonals(ns, nzrow, ds, offsets):
   """Build diagonals from function list and offsets"""

   diags = [0]*len(ds)
   for d,f in enumerate(ds):
      lb = max(0, -offsets[d])
      ub = min(len(ns),len(ns)-offsets[d])
      diags[d] = [f(n) if n > nzrow else 0 for n in ns[lb:ub]]

   return diags

