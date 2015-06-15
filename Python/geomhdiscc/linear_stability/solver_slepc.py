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

class GEVPSolver:
    """GEVP Solver using on SLEPc"""

    def __init__(self, shift_range = None, tol = 1e-8):
        """Initialize the SLEPc solver"""

        self.tol = tol

        self.create_eps(shift_range = None)

    def create_eps(self, shift_range = None):
        """Create SLEPc's eigensolver"""

        opts = PETSc.Options()
        opts["mat_mumps_icntl_14"] = 80
        opts["mat_mumps_icntl_29"] = 2

        if shift_range is None:
            #self.shift_range = (1e-2, 0.2)
            self.shift_range = (-1e-2, 1e-2)
        else:
            self.shift_range = shift_range

        self.E = SLEPc.EPS()
        self.E.create()

        self.E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
        #self.E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
        self.E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)
        self.E.setBalance(SLEPc.EPS.Balance.TWOSIDE)
        self.E.setTolerances(tol = self.tol)

        ST = self.E.getST()
        ST.setType('sinvert')
        if self.shift_range[0] == self.shift_range[1]:
            s = shift_range[0]
        else:
            rnd = PETSc.Random()
            rnd.create(comm = MPI.COMM_SELF)
            rnd.setType(PETSc.Random.Type.RAND)
            rnd.setInterval(self.shift_range)
            s = rnd.getValueReal()

        ST.setShift(s)
        KSP = ST.getKSP()
        KSP.setType('preonly')
        PC = KSP.getPC()
        PC.setType('lu')
        PC.setFactorSolverPackage('mumps')
        PC.setFromOptions()
        KSP.setFromOptions()
        ST.setFromOptions()

        self.E.setFromOptions()

    def update_eps(self, A, B, nev, initial_vector = None, shift_range = None):
        """Create SLEPc eigensolver"""

        self.create_eps(shift_range = shift_range)

        self.E.setOperators(A,B)
        if initial_vector is not None:
            v = PETSc.Vec().createWithArray(initial_vector)
            self.E.setInitialSpace(v)

        self.E.setDimensions(nev = nev)

#        ST = E.getST()
#        ST.setType('sinvert')
#        if shift_range[0] == shift_range[1]:
#            s = shift_range[0]
#        else:
#            rnd = PETSc.Random()
#            rnd.create(comm = MPI.COMM_SELF)
#            rnd.setType(PETSc.Random.Type.RAND)
#            rnd.setInterval(shift_range)
#            s = rnd.getValueReal()
#        ST.setShift(s)

    def eigenvalues(self, system, nev, initial_vector = None):
        """Compute eigenvalues using SLEPc"""

        pA, pB = self.petsc_operators(*system)

        self.update_eps(pA, pB, nev, initial_vector = initial_vector)

        self.E.solve()
        nconv = self.E.getConverged()

        if nconv >= nev:
            eigs = np.array(np.zeros(nev), dtype=complex)
            for i in range(nev):
                eigs[i] = self.E.getEigenvalue(i)

            return eigs
        else:
            return None

    def eigenpairs(self, system, nev, initial_vector = None):
        """Compute eigenpairs using SLEPc"""

        pA, pB = self.petsc_operators(*system)

        self.update_eps(pA, pB, nev, initial_vector = initial_vector)

        self.E.solve()
        nconv = self.E.getConverged()

        if nconv >= nev:
            # Create the results vectors
            vr, wr = pA.getVecs()
            vi, wi = pA.getVecs()
            eigs = np.array(np.zeros(nev), dtype=complex)
            vects = np.array(np.zeros((vi.getLocalSize(),nev)), dtype=complex)
            for i in range(nev):
                eigs[i] = self.E.getEigenpair(i, vr, vi)
                vects[:,i] = vr.getArray() + 1j*vi.getArray()

            return (eigs, vects)
        else:
            return (None, None)

    def restrict_operators(self, sizes):
        """Compute restriction for operators"""

        if MPI.COMM_WORLD.Get_size() > 1:
            pTmp = PETSc.Vec().create()
            pTmp.setSizes(np.sum(sizes[0]))
            pTmp.setUp()
            rstart, rend = pTmp.getOwnershipRange()
           
            restrict = []
            bstart = 0
            tot = 0
            for s, l in zip(sizes[0], sizes[1]):
                if rstart < tot + s:
                    bstart = max(np.floor((rstart-tot)/l),0) + sizes[2]
                    bend = min(np.ceil((rend-tot)/l),np.ceil(s/l)) + sizes[2]
                    restrict.append(np.arange(bstart,bend))
                else:
                    restrict.append(np.array([]))
                tot = tot + s

        else:
            restrict = None

        return restrict

    def petsc_operators(self, opA, opB, sizes):
        """Convert SciPy operators to PETSc operators"""

        # Build operator restriction
        #restrict = self.restrict_operators(sizes)
        restrict = None

        # Setup A matrix
        A = opA(restriction = restrict).transpose().tocsr()
        pA = PETSc.Mat().create()
        pA.setSizes(A.shape)
        pA.setUp()

        # Fill A matrix
        rstart, rend = pA.getOwnershipRange()
        pA.createAIJ(size=A.shape, nnz=A.getnnz(1)[rstart:rend+1], csr=(A.indptr[rstart:rend+1] - A.indptr[rstart], A.indices[A.indptr[rstart]:A.indptr[rend]], A.data[A.indptr[rstart]:A.indptr[rend]]))
        A = None
        pA.assemblyBegin()

        # Setup B matrix
        B = opB(restriction = restrict).transpose().tocsr()
        pA.assemblyEnd()
        pA.transpose()

        pB = PETSc.Mat().create()
        pB.setSizes(B.shape)
        pB.setUp()

        # Fill B matrix
        rstart, rend = pB.getOwnershipRange()
        pB.createAIJ(size=B.shape, nnz=B.getnnz(1)[rstart:rend+1], csr=(B.indptr[rstart:rend+1] - B.indptr[rstart], B.indices[B.indptr[rstart]:B.indptr[rend]], B.data[B.indptr[rstart]:B.indptr[rend]]))
        pB.assemble()
        B = None
        pB.transpose()

        return (pA, pB)

class NewtonNoneError(Exception):
    pass

class NewtonDoneError(Exception):
    pass

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
