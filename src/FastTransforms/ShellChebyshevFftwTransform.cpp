/** 
 * @file ShellChebyshevFftwTransform.cpp
 * @brief Source of the implementation of the Chebyshev FFTW transform for a spherical shell radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>
#include <limits>

// External includes
//

// Class include
//
#include "FastTransforms/ShellChebyshevFftwTransform.hpp"

// Configuration includes
//
#include "Profiler/ProfilerMacro.h"

// Project includes
//
#include "StaticAsserts/StaticAssert.hpp"
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"
#include "Python/PythonWrapper.hpp"
#include "FastTransforms/Validator/FftwTests.hpp"

namespace QuICC {

namespace Transform {

   Array ShellChebyshevFftwTransform::generateGrid(const int size, const MHDFloat ro, const MHDFloat rRatio)
   {
      if(ro > 0 && rRatio >= 0)
      {
         // Initialise grid storage
         Array grid(size);

         // Create Chebyshev grid
         for(int k = 0; k < size; k++)
         {
            grid(k) = std::cos((Math::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(size));

            MHDFloat b = (ro*rRatio + ro)/2.0;
            MHDFloat a = ro - b;

            grid(k) = a*grid(k) + b;
         }

         return grid;
      } else
      {
         throw Exception("generateGrid called with incompatible gap width or radii ratio");
      }
   }

   ShellChebyshevFftwTransform::ShellChebyshevFftwTransform()
      : mFPlan(NULL), mBPlan(NULL), mRo(-1), mRRatio(-1), mCnstA(0.0), mCnstB(0.0)
   {
      // Initialise the Python interpreter wrapper
      PythonWrapper::init();

      // Initialize FFTW
      FftwLibrary::initFft();
   }

   ShellChebyshevFftwTransform::~ShellChebyshevFftwTransform()
   {
      // Cleanup memory used by FFTW
      FftwLibrary::cleanupFft();
   }

   void ShellChebyshevFftwTransform::init(ShellChebyshevFftwTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Set the scaling factor
      this->mspSetup->setScale(1.0/static_cast<MHDFloat>(2*this->mspSetup->fwdSize()));

      // Initialise FFTW interface
      this->initFft();

      // Initialise Chebyshev operator(s)
      this->initOperators();

      // Register the FFTW object
      FftwLibrary::registerFft();
   }

   void ShellChebyshevFftwTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      list.insert(NonDimensional::RO);

      list.insert(NonDimensional::RRATIO);
   }

   void ShellChebyshevFftwTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      this->mRo = options.find(NonDimensional::RO)->second;

      this->mRRatio = options.find(NonDimensional::RRATIO)->second;

      // Compute the corresponding linear change of variable r = ax + b
      this->mCnstB = (this->mRo*this->mRRatio + this->mRo)/2.0;

      this->mCnstA = this->mRo - this->mCnstB;
   }

   Array ShellChebyshevFftwTransform::meshGrid() const
   {
      return ShellChebyshevFftwTransform::generateGrid(this->mspSetup->fwdSize(), this->mRo, this->mRRatio);
   }

   void ShellChebyshevFftwTransform::initFft()
   {  
      /// \mhdBug implement strided tranforms for complex <-> complex case if possible

      int fwdSize = this->mspSetup->fwdSize();
      int bwdSize = this->mspSetup->bwdSize();
      int howmany = this->mspSetup->howmany();

      // Create the two plans
      const int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::COMPONENT)
      {
      } else
      {
      }

      // Run validation tests
      Validator::FftwTests::validateR2R10();
      Validator::FftwTests::validateR2R01();

      // Initialise temporary storage
      this->mTmpIn.setZero(fwdSize, howmany);
      this->mTmpOut.setZero(bwdSize, howmany);

      // Create the physical to spectral plan
      const fftw_r2r_kind fwdKind[] = {FFTW_REDFT10};
      this->mFPlan = fftw_plan_many_r2r(1, fftSize, howmany, this->mTmpIn.data(), NULL, 1, fwdSize, this->mTmpOut.data(), NULL, 1, bwdSize, fwdKind, FftwLibrary::planFlag());

      // Create the spectral to physical plan
      const fftw_r2r_kind bwdKind[] = {FFTW_REDFT01};
      this->mBPlan = fftw_plan_many_r2r(1, fftSize, howmany, this->mTmpOut.data(), NULL, 1, bwdSize, this->mTmpIn.data(), NULL, 1, fwdSize, bwdKind, FftwLibrary::planFlag());
   }

   void ShellChebyshevFftwTransform::initOperators()
   {
      //
      // Initialise projector operators
      //

      //
      // Initialise integrator operators
      //
      // Multiplication by R
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGR, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // QST Q operator (4th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGQ4, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // QST S operator (4th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGS4, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // QST T operator
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGT, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // QST Q operator (2th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGQ2, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // QST S operator (2th order)
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGS2, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));

      //
      // Initialise solver operators
      //
      // Multiplication by R
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIVR, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // Multiplication by R^2
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIVR2, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // First derivative
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIFF, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // Second derivative
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIFF2, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));

      // Initialise python wrapper
      PythonWrapper::import("quicc.geometry.spherical.shell_radius");

      // Prepare arguments to Chebyshev matrices call
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(4);
      // ... compute a, b factors
      PyObject *pTmp = PyTuple_New(2);
      PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(this->mRo));
      PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(this->mRRatio));
      PythonWrapper::setFunction("linear_r2x");
      pValue = PythonWrapper::callFunction(pTmp);
      PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue, 0));
      PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue, 1));
      // ... create boundray condition (none)
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
      PyTuple_SetItem(pArgs, 3, pValue);

      // ... set forward operator size
      pValue = PyLong_FromLong(this->mspSetup->fwdSize());
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call r1
      PythonWrapper::setFunction("r1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGR)->second, pValue);
      Py_DECREF(pValue);

      // Call QST Q component (4th order)
      PythonWrapper::setFunction("i4r3");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGQ4)->second, pValue);
      Py_DECREF(pValue);

      // Call QST S component (4th order)
      PythonWrapper::setFunction("i4r3d1r1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGS4)->second, pValue);
      Py_DECREF(pValue);

      // Call QST T component
      PythonWrapper::setFunction("i2r2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGT)->second, pValue);
      Py_DECREF(pValue);

      // Call QST Q component (2nd order)
      PythonWrapper::setFunction("i2r1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGQ2)->second, pValue);
      Py_DECREF(pValue);

      // Call QST S component (2nd order)
      PythonWrapper::setFunction("i2r1d1r1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGS2)->second, pValue);
      Py_DECREF(pValue);

      // Change resoluton for solver matrices
      pValue = PyLong_FromLong(this->mspSetup->fwdSize());
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call r1 for solver
      PythonWrapper::setFunction("r1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mSolveOp.find(ProjectorType::DIVR)->second, pValue);
      Py_DECREF(pValue);
      // Rescale first row to create SPD matrix (due to c_0 in Chebyshev expansion)
      for(int j = 0; j < this->mSolveOp.find(ProjectorType::DIVR)->second.cols(); ++j)
      {
         if(this->mSolveOp.find(ProjectorType::DIVR)->second.coeff(0,j) != 0.0)
         {
            this->mSolveOp.find(ProjectorType::DIVR)->second.coeffRef(0,j) *= 0.5;
         }
      }

      // Call r2 for solver
      PythonWrapper::setFunction("r2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mSolveOp.find(ProjectorType::DIVR2)->second, pValue);
      Py_DECREF(pValue);
      // Rescale first row to create SPD matrix (due to c_0 in Chebyshev expansion)
      for(int j = 0; j < this->mSolveOp.find(ProjectorType::DIVR2)->second.cols(); ++j)
      {
         if(this->mSolveOp.find(ProjectorType::DIVR2)->second.coeff(0,j) != 0.0)
         {
            this->mSolveOp.find(ProjectorType::DIVR2)->second.coeffRef(0,j) *= 0.5;
         }
      }

      // Call i1 for solver
      // Change resolution
      pValue = PyLong_FromLong(this->mspSetup->fwdSize()+1);
      PyTuple_SetItem(pArgs, 0, pValue);
      // ... change boundary condition to zero last modes
      pValue = PyTuple_GetItem(pArgs, 3);
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
      PyDict_SetItem(pValue, PyUnicode_FromString("rt"), PyLong_FromLong(1));
      PyDict_SetItem(pValue, PyUnicode_FromString("cr"), PyLong_FromLong(1));
      PythonWrapper::setFunction("i1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mSolveOp.find(ProjectorType::DIFF)->second, pValue);
      Py_DECREF(pValue);

      // Call i2 for solver
      // Change resolution
      pValue = PyLong_FromLong(this->mspSetup->fwdSize()+2);
      PyTuple_SetItem(pArgs, 0, pValue);
      // ... change boundary condition to zero last modes
      pValue = PyTuple_GetItem(pArgs, 3);
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
      PyDict_SetItem(pValue, PyUnicode_FromString("rt"), PyLong_FromLong(2));
      PyDict_SetItem(pValue, PyUnicode_FromString("cr"), PyLong_FromLong(2));
      PythonWrapper::setFunction("i2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mSolveOp.find(ProjectorType::DIFF2)->second, pValue);
      Py_DECREF(pValue);

      // Initialize solver storage
      this->mTmpInS.setZero(this->mspSetup->fwdSize(), this->mspSetup->howmany());
      this->mTmpOutS.setZero(this->mspSetup->bwdSize(), this->mspSetup->howmany());

      // Initialize solver and factorize division by R operator (SPD)
      SharedPtrMacro<Solver::SparseSpdSelector<SparseMatrix>::Type>  pSpdSolver(new Solver::SparseSpdSelector<SparseMatrix>::Type());
      this->mSpdSolver.insert(std::make_pair(ProjectorType::DIVR, pSpdSolver));
      this->mSpdSolver.find(ProjectorType::DIVR)->second->compute(this->mSolveOp.find(ProjectorType::DIVR)->second);
      // Check for successful factorisation
      if(this->mSpdSolver.find(ProjectorType::DIVR)->second->info() != Eigen::Success)
      {
         throw Exception("Factorization of backward division by R failed!");
      }

      // Initialize solver and factorize division by R^2 operator (SPD)
      pSpdSolver = SharedPtrMacro<Solver::SparseSpdSelector<SparseMatrix>::Type>(new Solver::SparseSpdSelector<SparseMatrix>::Type());
      this->mSpdSolver.insert(std::make_pair(ProjectorType::DIVR2, pSpdSolver));
      this->mSpdSolver.find(ProjectorType::DIVR2)->second->compute(this->mSolveOp.find(ProjectorType::DIVR2)->second);
      // Check for successful factorisation
      if(this->mSpdSolver.find(ProjectorType::DIVR2)->second->info() != Eigen::Success)
      {
         throw Exception("Factorization of backward division by R^2 failed!");
      }

      // Initialize solver and factorize d^1 operator (upper triangular)
      SharedPtrMacro<Solver::SparseTriSelector<SparseMatrix>::Type>  pSolver(new Solver::SparseTriSelector<SparseMatrix>::Type());
      this->mTriSolver.insert(std::make_pair(ProjectorType::DIFF, pSolver));
      this->mTriSolver.find(ProjectorType::DIFF)->second->compute(this->mSolveOp.find(ProjectorType::DIFF)->second);
      // Check for successful factorisation
      if(this->mTriSolver.find(ProjectorType::DIFF)->second->info() != Eigen::Success)
      {
         throw Exception("Factorization of backward 1st derivative failed!");
      }

      // Initialize solver and factorize D^2 operator (upper triangular)
      pSolver = SharedPtrMacro<Solver::SparseTriSelector<SparseMatrix>::Type>(new Solver::SparseTriSelector<SparseMatrix>::Type());
      this->mTriSolver.insert(std::make_pair(ProjectorType::DIFF2, pSolver));
      this->mTriSolver.find(ProjectorType::DIFF2)->second->compute(this->mSolveOp.find(ProjectorType::DIFF2)->second);
      // Check for successful factorisation
      if(this->mTriSolver.find(ProjectorType::DIFF2)->second->info() != Eigen::Success)
      {
         throw Exception("Factorization of backward 2nd derivative failed!");
      }

      // Cleanup
      PythonWrapper::finalize();
   }

   void ShellChebyshevFftwTransform::cleanupFft()
   {
      // Detroy forward plan
      if(this->mFPlan)
      {
         fftw_destroy_plan(this->mFPlan);
      }

      // Detroy backward plan
      if(this->mBPlan)
      {
         fftw_destroy_plan(this->mBPlan);
      }

      // Unregister the FFTW object
      FftwLibrary::unregisterFft();

      // cleanup fftw library
      FftwLibrary::cleanupFft();
   }

   void ShellChebyshevFftwTransform::integrate(Matrix& rChebVal, const Matrix& physVal, ShellChebyshevFftwTransform::IntegratorType::Id integrator)
   {
      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Do transform
      DetailedProfilerMacro_start(ProfilerMacro::FWD1DTRAFFT);
      fftw_execute_r2r(this->mFPlan, const_cast<MHDFloat *>(physVal.data()), rChebVal.data());
      DetailedProfilerMacro_stop(ProfilerMacro::FWD1DTRAFFT);

      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTG)
      {
         // Rescale to remove FFT scaling
         rChebVal.topRows(this->mspSetup->specSize()) *= this->mspSetup->scale();

      } else
      {
         rChebVal.topRows(this->mspSetup->specSize()) = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second.topRows(this->mspSetup->specSize())*rChebVal;
      }

      #ifdef QUICC_DEBUG
         rChebVal.bottomRows(this->mspSetup->padSize()).setConstant(std::numeric_limits<MHDFloat>::quiet_NaN());
      #endif //QUICC_DEBUG
   }

   void ShellChebyshevFftwTransform::project(Matrix& rPhysVal, const Matrix& chebVal, ShellChebyshevFftwTransform::ProjectorType::Id projector)
   {
      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input matrix
      assert(chebVal.rows() == this->mspSetup->bwdSize());
      assert(chebVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-1) = chebVal.block(1, 0, this->mspSetup->specSize()-1, chebVal.cols()); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         this->mTmpIn = this->mTmpOutS;

      // Compute second derivative
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-2) = chebVal.block(2, 0, this->mspSetup->specSize()-2, chebVal.cols()); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+2).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF2);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF2);
         this->mTmpIn = this->mTmpOutS;

      // Compute division by R
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         // Rescale first row for SPD solve
         this->mTmpInS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSpdSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR);
         this->mTmpIn = this->mTmpOutS;

      // Compute division by R^2
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         // Rescale first row for SPD solve
         this->mTmpInS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR2);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSpdSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR2);
         this->mTmpIn = this->mTmpOutS;

      // Compute D r projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFFR)
      {
         this->mTmpInS.topRows(this->mTmpInS.rows()-1) = this->mSolveOp.find(ProjectorType::DIVR)->second.block(1, 0, this->mTmpInS.rows()-1, this->mspSetup->specSize())*chebVal.topRows(this->mspSetup->specSize());
         this->mTmpInS.bottomRows(1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         this->mTmpIn = this->mTmpOutS;

      // Compute 1/r D r projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         this->mTmpInS.topRows(this->mTmpInS.rows()-1) = this->mSolveOp.find(ProjectorType::DIVR)->second.block(1, 0, this->mTmpInS.rows()-1, this->mspSetup->specSize())*chebVal.topRows(this->mspSetup->specSize());
         this->mTmpInS.bottomRows(1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         // Rescale first row for SPD solve
         this->mTmpOutS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR);
         Solver::internal::solveWrapper(this->mTmpInS, *this->mSpdSolver.find(ProjectorType::DIVR)->second, this->mTmpOutS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR);
         this->mTmpIn = this->mTmpInS;

      // Compute radial laplacian projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::RADLAPL)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-1) = chebVal.block(1, 0, this->mspSetup->specSize()-1, chebVal.cols());
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         this->mTmpInS.topRows(this->mTmpInS.rows()-1) = this->mSolveOp.find(ProjectorType::DIVR2)->second.bottomRows(this->mTmpInS.rows()-1)*this->mTmpOutS;
         this->mTmpInS.bottomRows(1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         // Rescale first row for SPD solve
         this->mTmpOutS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR2);
         Solver::internal::solveWrapper(this->mTmpInS, *this->mSpdSolver.find(ProjectorType::DIVR2)->second, this->mTmpOutS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR2);
         this->mTmpIn = this->mTmpInS;

      // Compute simple projection
      } else
      {
         // Copy into other array
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize());
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
      }

      // Do transform
      DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAFFT);
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), rPhysVal.data());
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAFFT);
   }

   void ShellChebyshevFftwTransform::integrate(MatrixZ& rChebVal, const MatrixZ& physVal, ShellChebyshevFftwTransform::IntegratorType::Id integrator)
   {
      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // 
      // Transform real part
      //
      this->mTmpIn = physVal.real();

      DetailedProfilerMacro_start(ProfilerMacro::FWD1DTRAFFT);
      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());
      DetailedProfilerMacro_stop(ProfilerMacro::FWD1DTRAFFT);

      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTG)
      {
         rChebVal.topRows(this->mspSetup->specSize()).real() = this->mspSetup->scale()*this->mTmpOut.topRows(this->mspSetup->specSize());

      } else
      {
         rChebVal.topRows(this->mspSetup->specSize()).real() = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second.topRows(this->mspSetup->specSize())*this->mTmpOut;
      }

      // 
      // Transform imaginary part
      //
      this->mTmpIn = physVal.imag();

      DetailedProfilerMacro_start(ProfilerMacro::FWD1DTRAFFT);
      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());
      DetailedProfilerMacro_stop(ProfilerMacro::FWD1DTRAFFT);

      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTG)
      {
         rChebVal.topRows(this->mspSetup->specSize()).imag() = this->mspSetup->scale()*this->mTmpOut.topRows(this->mspSetup->specSize());

      } else
      {
         rChebVal.topRows(this->mspSetup->specSize()).imag() = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second.topRows(this->mspSetup->specSize())*this->mTmpOut;
      }

      #ifdef QUICC_DEBUG
         rChebVal.bottomRows(this->mspSetup->padSize()).setConstant(std::numeric_limits<MHDFloat>::quiet_NaN());
      #endif //QUICC_DEBUG
   }

   void ShellChebyshevFftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ShellChebyshevFftwTransform::ProjectorType::Id projector)
   {
      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input  matrix
      assert(chebVal.rows() == this->mspSetup->bwdSize());
      assert(chebVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative of real part
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-1) = chebVal.block(1, 0, this->mspSetup->specSize()-1, chebVal.cols()).real(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         this->mTmpIn = this->mTmpOutS;

      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-2) = chebVal.block(2, 0, this->mspSetup->specSize(), chebVal.cols()).real(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+2).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF2);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF2);
         this->mTmpIn = this->mTmpOutS;

      // Compute division by R of real part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).real(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         // Rescale first row for SPD solve
         this->mTmpInS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSpdSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR);
         this->mTmpIn = this->mTmpOutS;

      // Compute division by R^2 of real part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).real(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         // Rescale first row for SPD solve
         this->mTmpInS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR2);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSpdSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR2);
         this->mTmpIn = this->mTmpOutS;

      // Compute D r projection of real part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFFR)
      {
         this->mTmpInS.topRows(this->mTmpInS.rows()-1) = this->mSolveOp.find(ProjectorType::DIVR)->second.block(1, 0, this->mTmpInS.rows()-1, this->mspSetup->specSize())*chebVal.topRows(this->mspSetup->specSize()).real();
         this->mTmpInS.bottomRows(1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         this->mTmpIn = this->mTmpOutS;

      // Compute 1/r D r projection of real part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         this->mTmpInS.topRows(this->mTmpInS.rows()-1) = this->mSolveOp.find(ProjectorType::DIVR)->second.block(1, 0, this->mTmpInS.rows()-1, this->mspSetup->specSize())*chebVal.topRows(this->mspSetup->specSize()).real();
         this->mTmpInS.bottomRows(1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         // Rescale first row for SPD solve
         this->mTmpOutS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR);
         Solver::internal::solveWrapper(this->mTmpInS, *this->mSpdSolver.find(ProjectorType::DIVR)->second, this->mTmpOutS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR);
         this->mTmpIn = this->mTmpInS;

      // Compute radial laplacian projection of real part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::RADLAPL)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-1) = chebVal.block(1, 0, this->mspSetup->specSize()-1, chebVal.cols()).real();
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         this->mTmpInS.topRows(this->mTmpInS.rows()-1) = this->mSolveOp.find(ProjectorType::DIVR2)->second.bottomRows(this->mTmpInS.rows()-1)*this->mTmpOutS;
         this->mTmpInS.bottomRows(1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         // Rescale first row for SPD solve
         this->mTmpOutS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR2);
         Solver::internal::solveWrapper(this->mTmpInS, *this->mSpdSolver.find(ProjectorType::DIVR2)->second, this->mTmpOutS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR2);
         this->mTmpIn = this->mTmpInS;

      // Compute simple projection of real part
      } else
      {
         // Copy values into simple matrix
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).real();
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
      }

      // Do transform of real part
      DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAFFT);
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAFFT);
      rPhysVal.real() = this->mTmpOut;

      // Compute first derivative of imaginary part
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-1) = chebVal.block(1, 0, this->mspSetup->specSize()-1, chebVal.cols()).imag(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         this->mTmpIn = this->mTmpOutS;

      // Compute second derivative by R of imaginary part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-2) = chebVal.block(2, 0, this->mspSetup->specSize()-2, chebVal.cols()).imag(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+2).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF2);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF2);
         this->mTmpIn = this->mTmpOutS;

      // Compute division by R of imaginary part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).imag(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         // Rescale first row for SPD solve
         this->mTmpInS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSpdSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR);
         this->mTmpIn = this->mTmpOutS;

      // Compute division by R^2 of imaginary part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).imag(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         // Rescale first row for SPD solve
         this->mTmpInS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR2);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSpdSolver.find(projector)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR2);
         this->mTmpIn = this->mTmpOutS;

      // Compute D r projection of imaginary part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFFR)
      {
         this->mTmpInS.topRows(this->mTmpInS.rows()-1) = this->mSolveOp.find(ProjectorType::DIVR)->second.block(1, 0, this->mTmpInS.rows()-1, this->mspSetup->specSize())*chebVal.topRows(this->mspSetup->specSize()).imag();
         this->mTmpInS.bottomRows(1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         this->mTmpIn = this->mTmpOutS;

      // Compute 1/r D r projection of imaginary part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         this->mTmpInS.topRows(this->mTmpInS.rows()-1) = this->mSolveOp.find(ProjectorType::DIVR)->second.block(1, 0, this->mTmpInS.rows()-1, this->mspSetup->specSize())*chebVal.topRows(this->mspSetup->specSize()).imag();
         this->mTmpInS.bottomRows(1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         // Rescale first row for SPD solve
         this->mTmpOutS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR);
         Solver::internal::solveWrapper(this->mTmpInS, *this->mSpdSolver.find(ProjectorType::DIVR)->second, this->mTmpOutS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR);
         this->mTmpIn = this->mTmpInS;

      // Compute radial laplacian projection of imaginary part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::RADLAPL)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-1) = chebVal.block(1, 0, this->mspSetup->specSize()-1, chebVal.cols()).imag();
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         this->mTmpInS.topRows(this->mTmpInS.rows()-1) = this->mSolveOp.find(ProjectorType::DIVR2)->second.bottomRows(this->mTmpInS.rows()-1)*this->mTmpOutS;
         this->mTmpInS.bottomRows(1).setZero();
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRADIFF);
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(ProjectorType::DIFF)->second, this->mTmpInS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRADIFF);
         // Rescale first row for SPD solve
         this->mTmpOutS.topRows(1) *= 0.5;
         DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAR2);
         Solver::internal::solveWrapper(this->mTmpInS, *this->mSpdSolver.find(ProjectorType::DIVR2)->second, this->mTmpOutS);
         DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAR2);
         this->mTmpIn = this->mTmpInS;

      // Compute simple projection of imaginary part
      } else
      {
         // Rescale results
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).imag();
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
      }

      // Do transform of imaginary part
      DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRAFFT);
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRAFFT);
      rPhysVal.imag() = this->mTmpOut;
   }

   void ShellChebyshevFftwTransform::integrate_full(Matrix& rChebVal, const Matrix& physVal, ShellChebyshevFftwTransform::IntegratorType::Id integrator)
   {
      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Do transform
      fftw_execute_r2r(this->mFPlan, const_cast<MHDFloat *>(physVal.data()), rChebVal.data());

      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTG)
      {
         // Rescale to remove FFT scaling
         rChebVal *= this->mspSetup->scale();

      } else
      {
         rChebVal = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second*rChebVal;
      }
   }

   void ShellChebyshevFftwTransform::integrate_full(MatrixZ& rChebVal, const MatrixZ& physVal, ShellChebyshevFftwTransform::IntegratorType::Id integrator)
   {
      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      //
      // Transform real part
      //
      this->mTmpIn = physVal.real();

      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());

      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTG)
      {
         rChebVal.real() = this->mspSetup->scale()*this->mTmpOut;

      } else
      {
         rChebVal.real() = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second*this->mTmpOut;
      }

      //
      // Transform imaginary part
      //
      this->mTmpIn = physVal.imag();

      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());

      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTG)
      {
         rChebVal.imag() = this->mspSetup->scale()*this->mTmpOut;

      } else
      {
         rChebVal.imag() = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second*this->mTmpOut;
      }
   }

#ifdef QUICC_STORAGEPROFILE
   MHDFloat ShellChebyshevFftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}
