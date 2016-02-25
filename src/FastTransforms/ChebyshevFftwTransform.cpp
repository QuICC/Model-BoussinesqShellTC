/** 
 * @file ChebyshevFftwTransform.cpp
 * @brief Source of the implementation of the FFTW transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/ChebyshevFftwTransform.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array ChebyshevFftwTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create Chebyshev grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = std::cos((Math::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(size));
      }

      return grid;
   }

   ChebyshevFftwTransform::ChebyshevFftwTransform()
      : mFPlan(NULL), mBPlan(NULL), mCScale(0.0)
   {
      // Initialise the Python interpreter wrapper
      PythonWrapper::init();

      // Initialize FFTW
      FftwLibrary::initFft();
   }

   ChebyshevFftwTransform::~ChebyshevFftwTransform()
   {
      // Cleanup memory used by FFTW
      FftwLibrary::cleanupFft();
   }

   void ChebyshevFftwTransform::init(ChebyshevFftwTransform::SharedSetupType spSetup)
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

   void ChebyshevFftwTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      if(dimId == Dimensions::Transform::TRA1D)
      {
         list.insert(NonDimensional::SCALE1D);

      } else if(dimId == Dimensions::Transform::TRA2D)
      {
         list.insert(NonDimensional::SCALE2D);

      } else if(dimId == Dimensions::Transform::TRA3D)
      {
         list.insert(NonDimensional::SCALE3D);
      }
   }

   void ChebyshevFftwTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      if(dimId == Dimensions::Transform::TRA1D)
      {
         this->mCScale = options.find(NonDimensional::SCALE1D)->second;

      } else if(dimId == Dimensions::Transform::TRA2D)
      {
         this->mCScale = options.find(NonDimensional::SCALE2D)->second;

      } else if(dimId == Dimensions::Transform::TRA3D)
      {
         this->mCScale = options.find(NonDimensional::SCALE3D)->second;
      }
   }

   Array ChebyshevFftwTransform::meshGrid() const
   {
      return ChebyshevFftwTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void ChebyshevFftwTransform::initFft()
   {  
      /// \mhdBug implement strided stranforms for complex <-> complex case if possible

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

   void ChebyshevFftwTransform::initOperators()
   {
      //
      // Initialise projector operators
      //

      //
      // Initialise integrator operators
      //

      //
      // Initialise solver operators
      //
      //---------------------------------------------------------------
      // First derivative
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIFF, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));

      // Initialise python wrapper
      PythonWrapper::import("geomhdiscc.geometry.cartesian.cartesian_1d");

      // Prepare arguments to Chebyshev matrices call
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(3);
      // ... create boundray condition (last mode is zero)
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(991));
      PyTuple_SetItem(pArgs, 1, pValue);
      // ... set coefficient to 1.0
      pValue = PyFloat_FromDouble(1.0/this->mCScale);
      PyTuple_SetItem(pArgs, 2, pValue);

      // ... set resoluton for solver matrices
      pValue = PyLong_FromLong(this->mspSetup->fwdSize());
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call i1 for solver
      PythonWrapper::setFunction("i1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mSolveOp.find(ProjectorType::DIFF)->second, pValue);
      Py_DECREF(pValue);

      // Initialize solver storage
      this->mTmpInS.setZero(this->mspSetup->fwdSize(), this->mspSetup->howmany());
      this->mTmpOutS.setZero(this->mspSetup->bwdSize(), this->mspSetup->howmany());

      // Initialize solver and factorize division by d1 operator
      SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>  pSolver = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mSolver.insert(std::make_pair(ProjectorType::DIFF, pSolver));
      this->mSolver.find(ProjectorType::DIFF)->second->compute(this->mSolveOp.find(ProjectorType::DIFF)->second);
      // Check for successful factorisation
      if(this->mSolver.find(ProjectorType::DIFF)->second->info() != Eigen::Success)
      {
         throw Exception("Factorization of backward 1st derivative failed!");
      }

      //---------------------------------------------------------------
      // Second derivative
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIFF2, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));

      // Prepare arguments to Chebyshev matrices call
      pArgs = PyTuple_New(3);
      // ... create boundray condition (last mode is zero)
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(992));
      PyTuple_SetItem(pArgs, 1, pValue);
      // ... set coefficient to 1.0
      pValue = PyFloat_FromDouble(1.0/this->mCScale);
      PyTuple_SetItem(pArgs, 2, pValue);

      // ... set resoluton for solver matrices
      pValue = PyLong_FromLong(this->mspSetup->fwdSize());
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call i1 for solver
      PythonWrapper::setFunction("i2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mSolveOp.find(ProjectorType::DIFF2)->second, pValue);
      Py_DECREF(pValue);

      // Initialize solver storage
      this->mTmpInS.setZero(this->mspSetup->fwdSize(), this->mspSetup->howmany());
      this->mTmpOutS.setZero(this->mspSetup->bwdSize(), this->mspSetup->howmany());

      // Initialize solver and factorize division by d1 operator
      pSolver = SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>(new Solver::SparseSelector<SparseMatrix>::Type());
      this->mSolver.insert(std::make_pair(ProjectorType::DIFF2, pSolver));
      this->mSolver.find(ProjectorType::DIFF2)->second->compute(this->mSolveOp.find(ProjectorType::DIFF2)->second);
      // Check for successful factorisation
      if(this->mSolver.find(ProjectorType::DIFF2)->second->info() != Eigen::Success)
      {
         throw Exception("Factorization of backward 1st derivative failed!");
      }

      // Cleanup
      PythonWrapper::finalize();
   }

   void ChebyshevFftwTransform::cleanupFft()
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

   void ChebyshevFftwTransform::integrate(Matrix& rChebVal, const Matrix& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

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

      if(integrator == ChebyshevFftwTransform::IntegratorType::INTG)
      {
         // Rescale to remove FFT scaling
         rChebVal.topRows(this->mspSetup->specSize()) *= this->mspSetup->scale();

      } else
      {
         rChebVal.topRows(this->mspSetup->specSize()) = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second.topRows(this->mspSetup->specSize())*rChebVal;
      }

      #ifdef GEOMHDISCC_DEBUG
         rChebVal.bottomRows(this->mspSetup->padSize()).setConstant(std::numeric_limits<MHDFloat>::quiet_NaN());
      #endif //GEOMHDISCC_DEBUG
   }

   void ChebyshevFftwTransform::project(Matrix& rPhysVal, const Matrix& chebVal, ChebyshevFftwTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input  matrix
      assert(chebVal.rows() == this->mspSetup->bwdSize());
      assert(chebVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative
      if(projector == ChebyshevFftwTransform::ProjectorType::DIFF)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()); 
         this->mTmpInS.topRows(1).setZero();
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSolver.find(projector)->second, this->mTmpInS);
         this->mTmpIn = this->mTmpOutS;

//         // Recurrence relation
//         this->recurrenceDiff(this->mTmpIn, chebVal.topRows(this->mspSetup->specSize()));
//         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Compute second derivative
      } else if(projector == ChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()); 
         this->mTmpInS.topRows(2).setZero();
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSolver.find(projector)->second, this->mTmpInS);
         this->mTmpIn = this->mTmpOutS;

      // Compute simple projection
      } else
      {
         // Copy into other array
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize());
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
      }

      // Do transform
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), rPhysVal.data());
   }

   void ChebyshevFftwTransform::integrate(MatrixZ& rChebVal, const MatrixZ& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

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

      if(integrator == ChebyshevFftwTransform::IntegratorType::INTG)
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

      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());

      if(integrator == ChebyshevFftwTransform::IntegratorType::INTG)
      {
         rChebVal.topRows(this->mspSetup->specSize()).imag() = this->mspSetup->scale()*this->mTmpOut.topRows(this->mspSetup->specSize());

      } else
      {
         rChebVal.topRows(this->mspSetup->specSize()).imag() = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second.topRows(this->mspSetup->specSize())*this->mTmpOut;
      }

      #ifdef GEOMHDISCC_DEBUG
         rChebVal.bottomRows(this->mspSetup->padSize()).setConstant(std::numeric_limits<MHDFloat>::quiet_NaN());
      #endif //GEOMHDISCC_DEBUG
   }

   void ChebyshevFftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ChebyshevFftwTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

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
      if(projector == ChebyshevFftwTransform::ProjectorType::DIFF)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).real(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         this->mTmpInS.topRows(1).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSolver.find(projector)->second, this->mTmpInS);
         this->mTmpIn = this->mTmpOutS;

//         // Recurrence relation
//         this->recurrenceDiff(this->mTmpIn, chebVal.topRows(this->mspSetup->specSize()).real());
//         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Compute second derivative of real part
      } else if(projector == ChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).real(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         this->mTmpInS.topRows(2).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSolver.find(projector)->second, this->mTmpInS);
         this->mTmpIn = this->mTmpOutS;

      // Compute simple projection of real part
      } else
      {
         // Copy values into simple matrix
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).real();
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
      }

      // Do transform of real part
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
      rPhysVal.real() = this->mTmpOut;

      // Compute first derivative of imaginary part
      if(projector == ChebyshevFftwTransform::ProjectorType::DIFF)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).imag(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         this->mTmpInS.topRows(1).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSolver.find(projector)->second, this->mTmpInS);
         this->mTmpIn = this->mTmpOutS;

//         // Recurrence relation
//         this->recurrenceDiff(this->mTmpIn, chebVal.topRows(this->mspSetup->specSize()).imag());
//         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Compute second derivative of imaginary part
      } else if(projector == ChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).imag(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()).setZero();
         this->mTmpInS.topRows(2).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mSolver.find(projector)->second, this->mTmpInS);
         this->mTmpIn = this->mTmpOutS;

      // Compute simple projection of imaginary part
      } else
      {
         // Rescale results
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).imag();
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
      }

      // Do transform of imaginary part
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
      rPhysVal.imag() = this->mTmpOut;
   }

   void ChebyshevFftwTransform::integrate_full(Matrix& rChebVal, const Matrix& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

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

      if(integrator == ChebyshevFftwTransform::IntegratorType::INTG)
      {
         // Rescale to remove FFT scaling
         rChebVal *= this->mspSetup->scale();

      } else
      {
         rChebVal = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second*rChebVal;
      }
   }

   void ChebyshevFftwTransform::integrate_full(MatrixZ& rChebVal, const MatrixZ& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

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

      if(integrator == ChebyshevFftwTransform::IntegratorType::INTG)
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

      if(integrator == ChebyshevFftwTransform::IntegratorType::INTG)
      {
         rChebVal.imag() = this->mspSetup->scale()*this->mTmpOut;

      } else
      {
         rChebVal.imag() = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second*this->mTmpOut;
      }

   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat ChebyshevFftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
