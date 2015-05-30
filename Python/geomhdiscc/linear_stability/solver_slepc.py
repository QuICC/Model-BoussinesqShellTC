"""Module provides the functions to compute linear stability curves"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import numpy.linalg as nplin
import scipy.linalg as splin
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin

from petsc4py import PETSc
from slepc4py import SLEPc


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

def petsc_operators(A, B):
    """Convert SciPy operators to PETSc operators"""

    A = A.tocsr()
    B = B.tocsr()

    pA = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))
    pA.assemble()
    pB = PETSc.Mat().createAIJ(size=B.shape, csr=(B.indptr, B.indices, B.data))
    pB.assemble()

    return (pA, pB)

def slepc_eps(A, B, nev, sigma = 0.0):
    """Create SLEPc eigensolver"""

    opts = PETSc.Options()
    opts["mat_mumps_icntl_14"] = 80

    E = SLEPc.EPS()
    E.create()

    E.setOperators(A,B)
    E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
    E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
    #E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)
    E.setBalance(SLEPc.EPS.Balance.TWOSIDE)
    E.setDimensions(nev = nev, ncv = 50)
    ST = E.getST()
    ST.setType('sinvert')
    ST.setShift(sigma)
    KSP = ST.getKSP()
    E.setTolerances(max_it=200)
    PC = KSP.getPC()
    PC.setType('lu')
    PC.setFactorSolverPackage('mumps')

    E.setFromOptions()

    return E

def eigenvalues(A, B, nev):
    """Compute eigenvalues using SLEPc"""

    pA, pB = petsc_operators(A, B)

    E = slepc_eps(pA, pB, nev, 0.1)

    E.solve()
    nconv = E.getConverged()

    if nconv >= nev:
        eigs = np.array(np.zeros(nev), dtype=complex)
        err = np.array(np.zeros(nev), dtype=complex)
        for i in range(nev):
            eigs[i] = E.getEigenvalue(i)
            err[i] = E.computeRelativeError(i)

        return eigs

def eigenpairs(A, B, nev):
    """Compute eigenpairs using SLEPc"""

    pA, pB = petsc_operators(A, B)

    E = slepc_eps(pA, pB, nev, 0.1)

    E.solve()
    nconv = E.getConverged()

    if nconv >= nev:
        # Create the results vectors
        vr, wr = pA.getVecs()
        vi, wi = pA.getVecs()
        eigs = np.array(np.zeros(nev), dtype=complex)
        vects = np.array(np.zeros((A.shape[0],nev)), dtype=complex)
        err = np.array(np.zeros(nev), dtype=complex)
        for i in range(nev):
            eigs[i] = E.getEigenpair(i, vr, vi)
            vects[:,i] = vr.getArray() + 1j*vi.getArray()
            err[i] = E.computeRelativeError(i)

        return (eigs, vects)
