"""Module provides the functions to compute linear stability curves"""

from __future__ import division
from __future__ import unicode_literals

import random
import numpy as np
import numpy.linalg as nplin
import scipy.linalg as splin
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin

import sys, slepc4py
slepc4py.init(sys.argv)

from mpi4py import MPI

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

class NewtonNoneError(Exception):
    pass

class NewtonDoneError(Exception):
    pass

def petsc_operators(opA, opB):
    """Convert SciPy operators to PETSc operators"""

    import timeit
    start = timeit.default_timer()
#    if MPI.COMM_WORLD.Get_size() > 1:
#        # Make sure we have clean COO matrices
#        A = opA().tocoo()
#        A.sum_duplicates()
#
#        # Create PETSc A matrix
#        pA = PETSc.Mat().create()
#        pA.setSizes(A.shape)
#        pA.setUp()
#        # Fill A matrix
#        for i,j,v in zip(A.row, A.col, A.data):
#            pA.setValue(i,j,v)
#        A = None
#        pA.assemblyBegin()
#
#        B = opB().tocoo()
#        B.sum_duplicates()
#
#        pA.assemblyEnd()
#
#        # Setup PETSc B matrix
#        pB = PETSc.Mat().create()
#        pB.setSizes(B.shape)
#        pB.setUp()
#        # Fill B matrix
#        for i,j,v in zip(B.row, B.col, B.data):
#            pB.setValue(i,j,v)
#        pB.assemble()
#        B = None
#    else:
    # Setup A matrix
    A = opA().tocsr()
    pA = PETSc.Mat().create()
    pA.setSizes(A.shape)
    pA.setUp()

    # Fill A matrix
    rstart, rend = pA.getOwnershipRange()
    pA.createAIJ(size=A.shape, nnz=A.getnnz(1), csr=(A.indptr[rstart:rend+1] - A.indptr[rstart], A.indices[A.indptr[rstart]:A.indptr[rend]], A.data[A.indptr[rstart]:A.indptr[rend]]))
    A = None
    pA.assemble()

    # Setup B matrix
    B = opB().tocsr()
    pB = PETSc.Mat().create()
    pB.setSizes(B.shape)
    pB.setUp()

    # Fill B matrix
    rstart, rend = pB.getOwnershipRange()
    pB.createAIJ(size=B.shape, nnz=B.getnnz(1), csr=(B.indptr[rstart:rend+1] - B.indptr[rstart], B.indices[B.indptr[rstart]:B.indptr[rend]], B.data[B.indptr[rstart]:B.indptr[rend]]))
    pB.assemble()
    B = None

    Print("Operator construction time: {:g}".format(timeit.default_timer() - start))

    return (pA, pB)

def slepc_eps(A, B, nev, tracker = None, initial_vector = None):
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
    E.setDimensions(nev = nev)
    #E.setTolerances(max_it= 100)
    if initial_vector is not None:
        v = PETSc.Vec().createWithArray(initial_vector)
        E.setInitialSpace(v)

    ST = E.getST()
    ST.setType('sinvert')
    if tracker is None or abs(tracker) > 0.5:
        rnd = PETSc.Random()
        rnd.create(comm = MPI.COMM_SELF)
        rnd.setType(PETSc.Random.Type.RAND)
        rnd.setInterval((0.1, 0.5))
        s = rnd.getValueReal()
    else:
        s = 0.0
    ST.setShift(s)
    KSP = ST.getKSP()
    KSP.setType('preonly')
    PC = KSP.getPC()
    PC.setType('lu')
    PC.setFactorSolverPackage('mumps')
    PC.setFromOptions()
    KSP.setFromOptions()
    ST.setFromOptions()

    E.setFromOptions()

    return E

def eigenvalues(A, B, nev, tracker = None, initial_vector = None):
    """Compute eigenvalues using SLEPc"""

    pA, pB = petsc_operators(A, B)

    E = slepc_eps(pA, pB, nev, tracker, initial_vector = initial_vector)

    E.solve()
    nconv = E.getConverged()

    if nconv >= nev:
        eigs = np.array(np.zeros(nev), dtype=complex)
        err = np.array(np.zeros(nev), dtype=complex)
        for i in range(nev):
            eigs[i] = E.getEigenvalue(i)
            err[i] = E.computeRelativeError(i)

        return eigs
    else:
        return None

def eigenpairs(A, B, nev, tracker = None, initial_vector = None):
    """Compute eigenpairs using SLEPc"""

    pA, pB = petsc_operators(A, B)

    E = slepc_eps(pA, pB, nev, tracker, initial_vector = initial_vector)

    E.solve()
    nconv = E.getConverged()

    if nconv >= nev:
        # Create the results vectors
        vr, wr = pA.getVecs()
        vi, wi = pA.getVecs()
        eigs = np.array(np.zeros(nev), dtype=complex)
        vects = np.array(np.zeros((vi.getLocalSize(),nev)), dtype=complex)
        err = np.array(np.zeros(nev), dtype=complex)
        for i in range(nev):
            eigs[i] = E.getEigenpair(i, vr, vi)
            vects[:,i] = vr.getArray() + 1j*vi.getArray()
            err[i] = E.computeRelativeError(i)

        return (eigs, vects)
    else:
        return (None, None)

def newton(func, x0, args=(), tol=1.48e-8, maxiter=50, step = 1e-4):
    """
    Find a zero using secant method.
    Find a zero of the function `func` given a nearby starting point `x0`.
    Parameters
    ----------
    func : function
        The function whose zero is wanted. It must be a function of a
        single variable of the form f(x,a,b,c...), where a,b,c... are extra
        arguments that can be passed in the `args` parameter.
    x0 : float
        An initial estimate of the zero that should be somewhere near the
        actual zero.
    args : tuple, optional
        Extra arguments to be used in the function call.
    tol : float, optional
        The allowable error of the zero value.
    maxiter : int, optional
        Maximum number of iterations.
    Returns
    -------
    zero : float
        Estimated location where function is zero.
    """
    if tol <= 0:
        raise ValueError("tol too small (%g <= 0)" % tol)
    # Secant method
    p0 = x0
    q0 = func(*((p0,) + args))
    if q0 >= 0:
        p1 = x0*(1.0 - step) + step
    else:
        p1 = x0*(1.0 + step) + step
    q1 = func(*((p1,) + args))
    for iter in range(maxiter):
        if q1 == q0:
            if p1 != p0:
                msg = "Tolerance of %s reached" % (p1 - p0)
                warnings.warn(msg, RuntimeWarning)
            return (p1 + p0)/2.0
        else:
            p = p1 - q1*(p1 - p0)/(q1 - q0)
            if p < 0:
                factor = 15.0/10.0
                if q1 < 0:
                    p = p1*factor
                else:
                    p = p1/factor
        if abs(p - p1) < tol:
            return p
        p0 = p1
        q0 = q1
        p1 = p

        try:
            q1 = func(*((p1,) + args))
        except NewtonNoneError:
            if abs(p1/p0) > 10:
                p1 = 10*p0
                q1 = func(*((p1,) + args))
            elif abs(p1/p0) > 2:
                p1 = (p0 + p1)/2.0
                q1 = func(*((p1,) + args))
            else:
                raise
        except:
            raise
            
    msg = "Failed to converge after %d iterations, value is %s" % (maxiter, p)
    raise RuntimeError(msg)
