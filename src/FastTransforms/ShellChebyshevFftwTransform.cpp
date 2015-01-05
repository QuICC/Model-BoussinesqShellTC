/** 
 * @file ShellChebyshevFftwTransform.cpp
 * @brief Source of the implementation of the Chebyshev FFTW transform for a spherical shell radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/ShellChebyshevFftwTransform.hpp"

// Project includes
//
#include "StaticAsserts/StaticAssert.hpp"
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

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
      : mFPlan(NULL), mBPlan(NULL), mRo(-1), mRRatio(-1)
   {
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
      // First derivative
      this->mDiff.resize(this->mspSetup->specSize(),this->mspSetup->specSize());

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.spherical.shell_radius");

      #if defined GEOMHDISCC_TRANSOP_FORWARD
         // Initialise array for multiplication by R
         this->mR = this->meshGrid();
         // Initialise array for division by R
         this->mDivR = this->meshGrid().array().pow(-1);
         // Initialise array for division by R^2
         this->mDivR2 = this->meshGrid().array().pow(-2);

         // Prepare arguments to d1(...) call
         PyObject *pArgs, *pValue;
         pArgs = PyTuple_New(4);
         // ... get operator size
         pValue = PyLong_FromLong(this->mspSetup->specSize());
         PyTuple_SetItem(pArgs, 0, pValue);
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

         // Call d1
         PythonWrapper::setFunction("d1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->mDiff, pValue);
         Py_DECREF(pValue);

         // Call d2
         PythonWrapper::setFunction("d2");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->mDiff2, pValue);
         Py_DECREF(pValue);

      #elif defined GEOMHDISCC_TRANSOP_BACKWARD
         // Initialise matrix for division by R
         this->mDivR.resize(this->mspSetup->specSize(),this->mspSetup->specSize());
         // Initialise matrix for division by R
         this->mDivR2.resize(this->mspSetup->specSize(),this->mspSetup->specSize());

         // Prepare arguments to i1(...) call
         PyObject *pArgs, *pValue;
         pArgs = PyTuple_New(4);
         // ... get operator size
         pValue = PyLong_FromLong(this->mspSetup->specSize());
         PyTuple_SetItem(pArgs, 0, pValue);
         // ... compute a, b factors
         PyObject *pTmp = PyTuple_New(2);
         PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(this->mRo));
         PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(this->mRRatio));
         PythonWrapper::setFunction("linear_r2x");
         pValue = PythonWrapper::callFunction(pTmp);
         PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue, 0));
         PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue, 1));
         // ... create boundary condition (last mode zero)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(99));
         PyTuple_SetItem(pArgs, 3, pValue);

         // Call i1
         PythonWrapper::setFunction("i1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->mDiff, pValue);
         Py_DECREF(pValue);

         // Call i1
         PythonWrapper::setFunction("i2");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->mDiff2, pValue);
         Py_DECREF(pValue);

         // ... create boundary condition (none)
         pValue = PyTuple_GetItem(pArgs, 3);
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
         // Call x1
         PythonWrapper::setFunction("x1");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->mDivR, pValue);
         Py_DECREF(pValue);

         // Call x2
         PythonWrapper::setFunction("x2");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix
         PythonWrapper::fillMatrix(this->mDivR2, pValue);
         Py_DECREF(pValue);
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Cleanup
      PythonWrapper::finalize();

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
         // Factorize differentiation matrix and free memory
         this->mSDiff.compute(this->mDiff);
         // Check for successful factorisation
         if(this->mSDiff.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward first derivative failed!");
         }

         // Factorize differentiation matrix and free memory
         this->mSDiff2.compute(this->mDiff2);
         // Check for successful factorisation
         if(this->mSDiff2.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward first second failed!");
         }

         // Factorize division matrix and free memory
         this->mSDivR.compute(this->mDivR);
         // Check for successful factorisation
         if(this->mSDivR.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward division by R failed!");
         }

         // Factorize division matrix and free memory
         this->mSDivR2.compute(this->mDivR2);
         // Check for successful factorisation
         if(this->mSDivR2.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward division by R^2 failed!");
         }
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD
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

   void ShellChebyshevFftwTransform::integrate(Matrix& rChebVal, const Matrix& physVal, ShellChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
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

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute division by R
      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVR)
      {
         this->mTmpIn = this->mDivR.asDiagonal()*physVal;
      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGR)
      {
         this->mTmpIn = this->mR.asDiagonal()*physVal;
      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVRDIFFR)
      {
         this->mTmpIn = this->mDivR.asDiagonal()*physVal;
      } else
      {
         this->mTmpIn = physVal;
      }
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Do transform
      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), rChebVal.data());

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVRDIFFR)
      {
         this->mTmpIn = physVal;
         fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());
         rChebVal.topRows(this->mspSetup->specSize()) += this->mDiff*this->mTmpOut.topRows(this->mspSetup->specSize());
      }
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Rescale to remove FFT scaling
      rChebVal *= this->mspSetup->scale();

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute division by R
      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVR)
      {
         throw Exception("INTGDIVR operator is not yet implemented");
      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGR)
      {
         throw Exception("INTGDIVR operator is not yet implemented");
      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVRDIFFR)
      {
         throw Exception("INTGDIVRDIFFR operator is not yet implemented");
      }
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD
   }

   void ShellChebyshevFftwTransform::project(Matrix& rPhysVal, const Matrix& chebVal, ShellChebyshevFftwTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
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
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff*chebVal.topRows(this->mspSetup->specSize());
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiff, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Compute second derivative
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF2)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff2*chebVal.topRows(this->mspSetup->specSize());
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiff2, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute division by R
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR)
      {
         this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()); 
         Solver::internal::solveWrapper(this->mTmpOutS, this->mSDivR, this->mTmpInS);
         this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute division by R^2
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR2)
      {
         this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()); 
         Solver::internal::solveWrapper(this->mTmpOutS, this->mSDivR2, this->mTmpInS);
         this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute 1/r D r projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff*chebVal.topRows(this->mspSetup->specSize());

      #endif //defined GEOMHDISCC_TRANSOP_FORWARD
      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute radial laplacian projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::RADLAPL)
      {
         // Compute D^2 part
         this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff2*chebVal.topRows(this->mspSetup->specSize());

      #endif //defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute simple projection
      } else
      {
         // Copy into other array
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize());
      }

      // Set the padded values to zero
      this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), rPhysVal.data());

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute division by R
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR)
      {
         rPhysVal = this->mDivR.asDiagonal()*rPhysVal;

      // Compute 1/r^2 projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR2)
      {
         rPhysVal = this->mDivR2.asDiagonal()*rPhysVal;

      // Compute 1/r D r projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize());
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
         fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
         rPhysVal += this->mDivR.asDiagonal()*this->mTmpOut;

      // Compute laplacian projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::RADLAPL)
      {
         // Compute 2.0*D/r part
         this->mTmpIn.topRows(this->mspSetup->specSize()) = 2.0*this->mDiff*chebVal.topRows(this->mspSetup->specSize());
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
         fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
         rPhysVal += this->mDivR.asDiagonal()*this->mTmpOut;
      }
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute 1/r D r projection
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         throw Exception("DIVRDIFFR operator is not yet implemented");

      // Compute laplacian projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         throw Exception("LAPL operator is not yet implemented");
      }
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD
   }

   void ShellChebyshevFftwTransform::integrate(MatrixZ& rChebVal, const MatrixZ& physVal, ShellChebyshevFftwTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
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

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute integration of division by R of real part
      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVR)
      {
         this->mTmpIn = this->mDivR.asDiagonal()*physVal.real();

      // Compute integration of r of real part
      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGR)
      {
         this->mTmpIn = this->mR.asDiagonal()*physVal.real();

      // Compute integration of 1/r D r of real part
      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVRDIFFR)
      {
         this->mTmpIn = this->mDivR.asDiagonal()*physVal.real();

      // Compute transform of real part
      } else
      {
         this->mTmpIn = physVal.real();
      }
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());

      // Rescale FFT output
      rChebVal.real() = this->mspSetup->scale()*this->mTmpOut;

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVRDIFFR)
      {
         this->mTmpIn = physVal.real();
         fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());
         rChebVal.topRows(this->mspSetup->specSize()).real() += this->mspSetup->scale()*this->mDiff*this->mTmpOut.topRows(this->mspSetup->specSize());
      }
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute division by R
      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVR)
      {
         throw Exception("INTGDIVR operator is not yet implemented");

      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGR)
      {
         throw Exception("INTGR operator is not yet implemented");

      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVRDIFFR)
      {
         throw Exception("INTGDIVRDIFFR operator is not yet implemented");
      }
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute integration of division by R of imaginary part
      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVR)
      {
         this->mTmpIn = this->mDivR.asDiagonal()*physVal.imag();

      // Compute integration of r of imaginary part
      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGR)
      {
         this->mTmpIn = this->mR.asDiagonal()*physVal.imag();

      // Compute integration of 1/r D r of imaginary part
      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVRDIFFR)
      {
         this->mTmpIn = this->mDivR.asDiagonal()*physVal.imag();

      // Compute transform of imaginary part
      } else
      {
         this->mTmpIn = physVal.imag();
      }
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());

      // Rescale FFT output
      rChebVal.imag() = this->mspSetup->scale()*this->mTmpOut;

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVRDIFFR)
      {
         this->mTmpIn = physVal.imag();
         fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());
         rChebVal.topRows(this->mspSetup->specSize()).imag() += this->mspSetup->scale()*this->mDiff*this->mTmpOut.topRows(this->mspSetup->specSize());
      }
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute division by R
      if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVR)
      {
         throw Exception("INTGDIVR operator is not yet implemented");

      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGR)
      {
         throw Exception("INTGR operator is not yet implemented");

      } else if(integrator == ShellChebyshevFftwTransform::IntegratorType::INTGDIVRDIFFR)
      {
         throw Exception("INTGDIVRDIFFR operator is not yet implemented");
      }
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD
   }

   void ShellChebyshevFftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ShellChebyshevFftwTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
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
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).real();
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).real(); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiff, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute division by R of real part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR)
      {
         this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).real(); 
         Solver::internal::solveWrapper(this->mTmpOutS, this->mSDivR, this->mTmpInS);
         this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute division by R differentiation of real part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).real();
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Compute simple projection of real part
      } else
      {
         // Copy values into simple matrix
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).real();
      }

      // Set the padded values to zero
      this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform of real part
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
      rPhysVal.real() = this->mTmpOut;

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute division by R
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR)
      {
         rPhysVal.real() = this->mDivR.asDiagonal()*rPhysVal.real();

      // Compute 1/r D r projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).real();
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
         fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
         rPhysVal.real() += this->mDivR.asDiagonal()*this->mTmpOut;
      }
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute 1/r D r projection
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         throw Exception("DIVRDIFFR operator is not yet implemented");
      }
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD

      // Compute first derivative of imaginary part
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIFF)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).imag();
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).imag(); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiff, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute division by R of imaginary part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR)
      {
         this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).imag(); 
         Solver::internal::solveWrapper(this->mTmpOutS, this->mSDivR, this->mTmpInS);
         this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute division by R differentiation of imaginary part
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).imag();
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Compute simple projection of imaginary part
      } else
      {
         // Rescale results
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).imag();
      }

      // Set the padded values to zero
      this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform of imaginary part
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
      rPhysVal.imag() = this->mTmpOut;

      #if defined GEOMHDISCC_TRANSOP_FORWARD
      // Compute division by R
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVR)
      {
         rPhysVal.imag() = this->mDivR.asDiagonal()*rPhysVal.imag();

      // Compute 1/r D r projection
      } else if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).imag();
         this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();
         fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
         rPhysVal.imag() += this->mDivR.asDiagonal()*this->mTmpOut;
      }
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
      // Compute 1/r D r projection
      if(projector == ShellChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
      {
         throw Exception("DIVRDIFFR operator is not yet implemented");
      }
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat ShellChebyshevFftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
