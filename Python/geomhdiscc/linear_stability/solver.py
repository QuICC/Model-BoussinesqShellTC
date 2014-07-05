"""Module provides the functions to compute linear stability curves"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import numpy.linalg as nplin
import scipy.linalg as splin
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin


def bound_gevp(A, B, mode, lb, ub):
    """Compute eigenvalues and bounds"""

    # Compute eigenvalues
    evp_vec, evp_lmb = solve_gevp(A, B, lb, ub)

    # Found an upper bound for eigenvalue search
    if sum(np.real(evp_lmb) > 0) >= mode:
        bound = evp_lmb(mode-1)
        lmb = evp_lmb(mode-1)
    elif len(evp_lmb) >= mode:
        bound = np.real(evp_lmb(min(5,len(evp_lmb))))
        lmb = evp_lmb(mode)
    else:
        bound = -1e0
        lmb = -np.inf

    return (bound, lmb)


def solve_gevp(A, B, lb, ub):
    """Solve the generalized eigenvalue problem"""

    # max multiplicity of eigenvalues
    max_mult = 100
    # restart whole calculation if it fails (random algorithm)
    max_fail = 10

    # basis size
    nbasis = 200

    # Loop with increasing multiplicity until all eigenvalue in interval have been found
    iresult = -1
    for mult in range(10, min(max_mult, A.shape[0]), 10):
        for run in range(0, max_fail):
            try:
                evp_vec, evp_lmb, iresult = sptarn(A, B, lb, ub, 100*np.spacing(2), nbasis, mult)
                break
            except:
                print("SPTARN crashed, restarting with different guess")

        if iresult > -1:
            break

    # Sort the eigenvalues and remove infinities
    evp_vec, evp_lmb = sort_no_inf(evp_vec, evp_lmb)

    return (evp_vec, evp_lmb)


def sort_no_inf(vec, lmb):
    """Sort eigenvalues and eigenvectors and remove infinities"""

    idx = np.argsort(np.real(lmb))[::-1]

    return (vec[:,idx], lmb[idx])


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
