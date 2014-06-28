"""Module provides the functions to compute linear stability curves"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import numpy.linalg as nplin
import scipy.linalg as splin
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import os.path

def trace_marginal(pbm, res, ks, params, save_file = True, new_file = False, db_file = False):
   """Compute the marginal curve for a given problem"""

   # Check for file output
   if save_file:
      # Exact filename is provided
      if isinstance(save_file, str):
         filename = save_file
         save_file = True
      # Use default naming scheme
      else:
         filename = 'marginal_' + pbm + '.dat'
      
      # Open file
      if new_file:
         f = open(filename, 'w')
      else:
         f = open(filename, 'a')
   else:
      points = []

   # Write file header if required
   if save_file and (new_file or os.path.isfile(filename)):
      write_header(f, pbm, res, params)

   # Compute a unique index
   tot = 1
   key_size = {}
   for k,v in params.items():
      try:
         sze = len(v)
      except:
         sze = 1
      tot = sze*tot
      key_size[k] = sze

   # Loop over all parameters
   for i in range(tot):
      # Create parameters for current run
      idx = dict(zip(list(key_size.keys()),np.unravel_index(i, tuple(key_size.values()))))
      run = {}
      for k,v in params.items():
         try:
            val = v[idx[k]]
         except:
            val = v
         run[k] = val

      # Loop over the wave numbers
      for kc in ks:
         # Compute critical values
         racs, omegas = find_critical(pbm, res, kc, rac0, run) 

         # Print run information
         if False:
            print('\t--> Converged to critical value: ')
            print('\t--> Converged to critical frequency: ')
         else:
            print('\t--> Failed to converge to critical value!')


         # Save results to file is required
         if save_file:
            write_results(f, res, kc, racs, omegas, run)
         else:
            if True:
               points.append([kc]+ racs + omegas)


   # Close file if required
   if save_file:
      f.close()
   else:
      return points


def find_critical(pbm, res, kc, rac0, params):
   """Find the critical value for given parameters"""

   # Loop over the requested modes
   for n in range(1,params['mode']+1):
      # Compute approximation with Brent's root finding algorithm
      print('Start root finding algorithm')


def write_header(f, pbm, res, params):
   """Write marginal curve file header"""
   
   # First header
   f.write('#Results for: ' + pbm + '\n')

   header = []

   # Resolution
   for i,r in enumerate(res):
      header.append('Res_' + str(i)) 

   # Parameters
   for k,v in sorted(params.items()):
      header.append(k)

   # Wave number
   header.append('k')

   # Critical values
   for i in range(1, params['mode']+1):
      header.append('Rac_' + str(i))

   # Critical frequencies
   for i in range(1, params['mode']+1):
      header.append('Omega_' + str(i))

   f.write('#'+'\t'.join(header) + '\n')


def write_results(f, res, kc, racs, omegas, params):
   """Write marginal curve point to file"""

   result = []

   # Resolution
   for r in res:
      result.append(str(r)) 

   # Parameters
   for k,v in sorted(params.items()):
      result.append(str(v))

   # Wave number
   result.append(str(kc))

   # Critical values
   for rac in racs:
      result.append(str(rac))

   # Critical frequencies
   for omega in omegas:
      result.append(str(omega))

   f.write('\t'.join(result) + '\n')


def sptarn(A, B, lb, ub, tolconv = 100*np.spacing(1), jmax = 100, maxmul = 10):
   """Compute eigenvalues in a given interval for the generalized eigenvalue problem using the shift-invert Arnoldi algorithm"""

   mode = 0

   # Convert sparse matrix to suitable format
   A = A.tocsc()
   B = B.tocsc()

   # Get matrix size
   n = A.shape[0]

   # Choose initial shift
   shift = choose_shift(lb, ub)

   # Limit number of basis vector to matrix size
   if jmax > n:
      jmax = n

   # Factorize operator and check for singularity
   maxtries = 3
   ntries = 0
   while ntries < maxtries:
      # Build shift-invert operator and factorize
      C = A - shift*B
      F = spsplin.factorized(C)

      # Check for singularity
      t = B*np.random.rand(n,1)
      t = t/nplin.norm(t,np.inf)
      v = F(t)
      fact = nplin.norm(v,np.inf)*np.max(np.abs(C)*np.ones((n,1)))

      # If factor is too large try new shift
      if fact > 1e-2/np.spacing(1):
         shift = choose_shift(lb, ub, True)

      else:
         del C
         break

      ntries = ntries + 1

   print('STARTING')
   # Could not find a suitable shift
   if ntries == maxtries:
      raise RuntimeError("Could not find a suitable shift!")

   # Number of eigenvalues with positive real part
   nPos = 0
   # Tolerances
   tol = tolconv
   tolspur = 100*tolconv
   tolstab = 1e-5
   # Convergence check points: after jf steps, then every js steps
   js = np.ceil(4000/n)
   jf = 10

   normhj = 0
   # nt number of converged values (inside and outside)
   nt = 0
   # US is n*nt matrix of Schur vectors
   US = np.empty((n,0))

   # Starting vector
   vstart = F(B*np.random.randn(n,1))
   vstart = F(B*np.ones((n,1)))

   terminatealg = False
   ################## CLEAN UP TO HERE ######################
   for mul in range(maxmul):
      # Orthogonalize starting vector against converged vectors
      v = np.hstack((US, ortho_gs(vstart, US)[0]))

      try:
         h = np.vstack((T, np.zeros((1,nt))))
      except NameError:
         h = np.zeros((1,nt))

      convrun = False

      # First testing point
      jt = min(nt+jf, jmax)
      for j in range(nt+1,jmax+1):
         r = F(B*(v[:,j-1:j]))
         vj1, hj = ortho_gs(r, v)
         h = np.hstack((h, hj[0:j]))
         normhj = max(normhj, nplin.norm(hj))
         # Test for linear dependence
         if hj[j] < tol*normhj:
            convrun = True
            jt = j
         else:
            v = np.hstack((v, vj1))
         
         h = np.vstack((h, np.hstack((np.zeros((1,j-1)), hj[j:j+1,:]))))

         # One vector added, time to test?
         if j == jt:
            iv = np.s_[nt:jt]
            ts, us = splin.schur(h[iv,iv], output='real')
            tc = splin.rsf2csf(ts, us)[0]
            lms = shift + 1/np.diag(tc)
            sij = np.abs((h[j,j-1]/normhj)*us[j-nt-1,:]).T

            convoutside = np.logical_and(np.logical_or(lms <= lb, lms > ub) , sij <= tolstab)
            nonconvinside = np.logical_and(np.logical_and(lb < lms, lms <= ub), np.logical_and(sij > tol, np.abs(us[0,:]).T > tolspur))
            convrun = convrun or (np.any(convoutside) and not np.any(nonconvinside))

            if np.any(nonconvinside):
               incinside = np.nonzero(nonconvinside)[0]
               iintr = incinside[0]
            else:
               iintr = 1

            print('\t\tBasis = ' + str(j) + ' New converged eig = ' + str(np.sum(sij<=tol)))
            jt = min(jt+js,jmax)
         if convrun:
            break

      nc = np.sum(np.cumprod(sij<=tol))
      print('End of sweep:\tBasis = ' + str(j) + ' New converged eig = ' + str(nc))

      if nc > 0:
         try:
            T = np.vstack((np.hstack((T, h[0:nt,iv].dot(us[:,0:nc]))), np.hstack((np.zeros((nc,nt)), ts[0:nc,0:nc]))))
         except NameError:
            T = np.vstack((h[0:nt,iv].dot(us[:,0:nc]), np.hstack((np.zeros((nc,nt)), ts[0:nc,0:nc]))))
         US = np.hstack((US, v[:,iv].dot(us[:,0:nc])))
         nt = nt + nc
         if nt >= jmax:
            terminatealg = True


      newinside = np.logical_and(lb < lms ,lms <= lb)
      if nc > 0 and np.any(np.logical_and(newinside ,sij <= tol)):
         vstart = F(B*np.random.randn(n,1))
      elif np.any(newinside):
         vstart = v[:,iv].dot(us[:,np.nonzero(newinside)[0][0]])
      else:
         terminatealg = True

      if ub == np.inf and np.any(np.logical_and(np.logical_and(newinside ,sij <= tol), np.real(lms) > 0)):
         nPos = nPos + np.sum(np.logical_and(np.logical_and(newinside, sij < tol) , np.real(lms) > 0))
         if mode > 0 and nPos > mode:
            terminatealg = True
            print('Early termination: Not all eigenvalues have been computed!')
      
      if terminatealg:
         break
   
   dt, s = nplin.eig(np.multiply(T,(np.abs(T)>tol*normhj)))
   ev = shift + 1.0/dt

   inside = np.nonzero(np.logical_and(lb < np.real(ev) ,np.real(ev) <= ub))[0]
   iresult = len(inside)
   if iresult > 0:
      ev = ev[inside]
      ievsrt = np.argsort(ev)
      lmb = ev[ievsrt]
      xv = US.dot(s[:,inside[ievsrt]])
   else:
      lmb = np.array([])
      xv = np.array([[]])

   if not terminatealg:
      iresult = -iresult

   return (xv, lmb, iresult)


def choose_shift(lb, ub, perturb = False):
   """Choose a shift in the [lb, ub] interval"""

   # No information about eigenvalue location
   if lb == -np.inf and ub == np.inf:
      shift = np.random.randn()

   # Use upper bound as shift
   elif lb == -np.inf:
      if perturb:
         # No information about eigenvalue location
         if ub == 0:
            shift = -np.abs(np.random.randn())
         else:
            shift = ub - abs(lb)*np.random.rand()

      else:
         shift = ub

   # Use lower bound as shift
   elif ub == np.inf:
      if perturb:
         # No information about eigenvalue location
         if lb == 0:
            shift = np.abs(np.random.randn())
         else:
            shift = lb + np.abs(lb)*np.random.rand()
      else:
         shift = lb

   # Pick random shift within range
   else:
      shift = (ub - lb)*np.random.rand() + lb

   return shift


def ortho_gs(v, mat):
   """Gram-Schmidt orthogonalization of vector against columns of matrix"""

   # Maximum number of reorthogonalizations
   maxort = 4
   # Orthogonaliation quality factor
   qmax = 0.7

   n, j = mat.shape
   vnn = nplin.norm(v)
   hc = np.zeros((j,1))

   if j > 0:
      for iort in range(maxort):
         hv = mat.conj().T.dot(v)
         hc = hc[:,0:1] + hv[:,0:1]
         v = v - mat.dot(hv)
         vnorm = vnn
         vnn = nplin.norm(v)
         # Check for reorthogonalization
         if vnn > qmax*vnorm:
            break

   # Build orthogonalization coefficients
   h = np.vstack((hc, vnn))
   y = np.reshape(v/vnn, (-1,1))

   return (y, h)
