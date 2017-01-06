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
      // i2 integrator
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI2, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // i4 integrator
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI4, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // i2d1 integrator
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI2D1, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // i4d1 integrator
      this->mIntgOp.insert(std::make_pair(IntegratorType::INTGI4D1, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));

      //
      // Initialise solver operators
      //
      //---------------------------------------------------------------
      // First derivative
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIFF, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));
      // Second derivative
      this->mSolveOp.insert(std::make_pair(ProjectorType::DIFF2, SparseMatrix(this->mspSetup->fwdSize(),this->mspSetup->fwdSize())));

      // Initialise python wrapper
      PythonWrapper::import("quicc.geometry.cartesian.cartesian_1d");

      // Prepare arguments to Chebyshev matrices call
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(3);
      // ... create boundray condition (none)
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
      PyTuple_SetItem(pArgs, 1, pValue);
      // ... set coefficient to 1.0
      pValue = PyFloat_FromDouble(1.0);
      PyTuple_SetItem(pArgs, 2, pValue);

      // ... set resolution for solver matrices
      pValue = PyLong_FromLong(this->mspSetup->fwdSize());
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call i2
      PythonWrapper::setFunction("i2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGI2)->second, pValue);
      Py_DECREF(pValue);

      // Call i4
      PythonWrapper::setFunction("i4");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGI4)->second, pValue);
      Py_DECREF(pValue);

      // ... set coefficient to cscale
      pValue = PyFloat_FromDouble(this->mCScale);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Call i2d1
      PythonWrapper::setFunction("i2d1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGI2D1)->second, pValue);
      Py_DECREF(pValue);

      // Call i4d1
      PythonWrapper::setFunction("i4d1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix
      PythonWrapper::fillMatrix(this->mIntgOp.find(IntegratorType::INTGI4D1)->second, pValue);
      Py_DECREF(pValue);

      // Call i1 for solver
      // Change resolution
      pValue = PyLong_FromLong(this->mspSetup->fwdSize()+1);
      PyTuple_SetItem(pArgs, 0, pValue);
      // ... set coefficient to cscale
      pValue = PyFloat_FromDouble(1.0/this->mCScale);
      PyTuple_SetItem(pArgs, 2, pValue);
      // ... change boundary condition to zero last modes
      pValue = PyTuple_GetItem(pArgs, 1);
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
      // ... set coefficient to cscale
      pValue = PyFloat_FromDouble(1.0/(this->mCScale*this->mCScale));
      PyTuple_SetItem(pArgs, 2, pValue);
      // ... change boundary condition to zero last modes
      pValue = PyTuple_GetItem(pArgs, 1);
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

      // Initialize solver and factorize d^1 operator (upper triangular)
      SharedPtrMacro<Solver::SparseTriSelector<SparseMatrix>::Type>  pSolver = SharedPtrMacro<Solver::SparseTriSelector<SparseMatrix>::Type>(new Solver::SparseTriSelector<SparseMatrix>::Type());
      this->mTriSolver.insert(std::make_pair(ProjectorType::DIFF, pSolver));
      this->mTriSolver.find(ProjectorType::DIFF)->second->compute(this->mSolveOp.find(ProjectorType::DIFF)->second);
      // Check for successful factorisation
      if(this->mTriSolver.find(ProjectorType::DIFF)->second->info() != Eigen::Success)
      {
         throw Exception("Factorization of backward 1st derivative failed!");
      }

      // Initialize solver and factorize d^2 operator (upper triangular)
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

   void ChebyshevFftwTransform::integrate(Matrix& rChebVal, const Matrix& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator)
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

      if(integrator == ChebyshevFftwTransform::IntegratorType::INTG)
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

   void ChebyshevFftwTransform::project(Matrix& rPhysVal, const Matrix& chebVal, ChebyshevFftwTransform::ProjectorType::Id projector)
   {
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
         this->mTmpInS.topRows(this->mspSetup->specSize()-1) = chebVal.block(1, 0, this->mspSetup->specSize()-1, chebVal.cols()); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+1).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
         this->mTmpIn = this->mTmpOutS;

      // Compute second derivative
      } else if(projector == ChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-2) = chebVal.block(2, 0, this->mspSetup->specSize()-2, chebVal.cols()); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+2).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
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

   void ChebyshevFftwTransform::integrate(MatrixZ& rChebVal, const MatrixZ& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator)
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

      if(integrator == ChebyshevFftwTransform::IntegratorType::INTG)
      {
         rChebVal.topRows(this->mspSetup->specSize()).real() = this->mspSetup->scale()*this->mTmpOut.topRows(this->mspSetup->specSize());

      } else if(integrator == ChebyshevFftwTransform::IntegratorType::INTGI2D1MI2)
      {
         // Mean has different operator
         if(this->mspSetup->idBlocks().size() > 0)
         {
            rChebVal.block(0, 0, this->mspSetup->specSize(), 1).real() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI2)->second.topRows(this->mspSetup->specSize())*this->mTmpOut.col(0);
            rChebVal.block(0, 1, this->mspSetup->specSize(), rChebVal.cols()-1).real() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI2D1)->second.topRows(this->mspSetup->specSize())*this->mTmpOut.rightCols(rChebVal.cols()-1);
         } else
         {
            rChebVal.topRows(this->mspSetup->specSize()).real() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI2D1)->second.topRows(this->mspSetup->specSize())*this->mTmpOut;
         }

      } else if(integrator == ChebyshevFftwTransform::IntegratorType::INTGI4D1MI2)
      {
         // Mean has different operator
         if(this->mspSetup->idBlocks().size() > 0)
         {
            rChebVal.block(0, 0, this->mspSetup->specSize(), 1).real() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI2)->second.topRows(this->mspSetup->specSize())*this->mTmpOut.col(0);
            rChebVal.block(0, 1, this->mspSetup->specSize(), rChebVal.cols()-1).real() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI4D1)->second.topRows(this->mspSetup->specSize())*this->mTmpOut.rightCols(rChebVal.cols()-1);
         } else
         {
            rChebVal.topRows(this->mspSetup->specSize()).real() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI4D1)->second.topRows(this->mspSetup->specSize())*this->mTmpOut;
         }

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

      } else if(integrator == ChebyshevFftwTransform::IntegratorType::INTGI2D1MI2)
      {
         // Mean has different operator
         if(this->mspSetup->idBlocks().size() > 0)
         {
            rChebVal.block(0, 0, this->mspSetup->specSize(), 1).imag() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI2)->second.topRows(this->mspSetup->specSize())*this->mTmpOut.col(0);
            rChebVal.block(0, 1, this->mspSetup->specSize(), rChebVal.cols()-1).imag() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI2D1)->second.topRows(this->mspSetup->specSize())*this->mTmpOut.rightCols(rChebVal.cols()-1);
         } else
         {
            rChebVal.topRows(this->mspSetup->specSize()).imag() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI2D1)->second.topRows(this->mspSetup->specSize())*this->mTmpOut;
         }

      } else if(integrator == ChebyshevFftwTransform::IntegratorType::INTGI4D1MI2)
      {
         // Mean has different operator
         if(this->mspSetup->idBlocks().size() > 0)
         {
            rChebVal.block(0, 0, this->mspSetup->specSize(), 1).imag() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI2)->second.topRows(this->mspSetup->specSize())*this->mTmpOut.col(0);
            rChebVal.block(0, 1, this->mspSetup->specSize(), rChebVal.cols()-1).imag() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI4D1)->second.topRows(this->mspSetup->specSize())*this->mTmpOut.rightCols(rChebVal.cols()-1);
         } else
         {
            rChebVal.topRows(this->mspSetup->specSize()).imag() = this->mspSetup->scale()*this->mIntgOp.find(ChebyshevFftwTransform::IntegratorType::INTGI4D1)->second.topRows(this->mspSetup->specSize())*this->mTmpOut;
         }

      } else
      {
         rChebVal.topRows(this->mspSetup->specSize()).imag() = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second.topRows(this->mspSetup->specSize())*this->mTmpOut;
      }

      #ifdef QUICC_DEBUG
         rChebVal.bottomRows(this->mspSetup->padSize()).setConstant(std::numeric_limits<MHDFloat>::quiet_NaN());
      #endif //QUICC_DEBUG
   }

   void ChebyshevFftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ChebyshevFftwTransform::ProjectorType::Id projector)
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
      if(projector == ChebyshevFftwTransform::ProjectorType::DIFF)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-1) = chebVal.block(1, 0, this->mspSetup->specSize()-1, chebVal.cols()).real(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+1).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
         this->mTmpIn = this->mTmpOutS;

      // Compute second derivative of real part
      } else if(projector == ChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-2) = chebVal.block(2, 0, this->mspSetup->specSize()-2, chebVal.cols()).real(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+2).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
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
         this->mTmpInS.topRows(this->mspSetup->specSize()-1) = chebVal.block(1, 0, this->mspSetup->specSize()-1, chebVal.cols()).imag(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+1).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
         this->mTmpIn = this->mTmpOutS;

      // Compute second derivative of imaginary part
      } else if(projector == ChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         this->mTmpInS.topRows(this->mspSetup->specSize()-2) = chebVal.block(2, 0, this->mspSetup->specSize()-2, chebVal.cols()).imag(); 
         this->mTmpInS.bottomRows(this->mspSetup->padSize()+2).setZero();
         Solver::internal::solveWrapper(this->mTmpOutS, *this->mTriSolver.find(projector)->second, this->mTmpInS);
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

   void ChebyshevFftwTransform::integrate_full(Matrix& rChebVal, const Matrix& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator)
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

      if(integrator == ChebyshevFftwTransform::IntegratorType::INTG)
      {
         // Rescale to remove FFT scaling
         rChebVal *= this->mspSetup->scale();

      } else
      {
         rChebVal = this->mspSetup->scale()*this->mIntgOp.find(integrator)->second*rChebVal;
      }
   }

   void ChebyshevFftwTransform::integrate_full(MatrixZ& rChebVal, const MatrixZ& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator)
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

#ifdef QUICC_STORAGEPROFILE
   MHDFloat ChebyshevFftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}
