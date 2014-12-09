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

      #if defined GEOMHDISCC_TRANSOP_FORWARD || defined GEOMHDISCC_TRANSOP_BACKWARD
         // Initialise Chebyshev operator(s)
         this->initOperators();
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD || defined GEOMHDISCC_TRANSOP_BACKWARD

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

   #if defined GEOMHDISCC_TRANSOP_FORWARD || defined GEOMHDISCC_TRANSOP_BACKWARD
   void ChebyshevFftwTransform::initOperators()
   {
      // Storage for the differentiation operator
      this->mDiff.resize(this->mspSetup->specSize(),this->mspSetup->specSize());

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.cartesian.cartesian_1d");

      #if defined GEOMHDISCC_TRANSOP_FORWARD
         // Prepare arguments to d1(...) call
         PyObject *pArgs, *pValue;
         pArgs = PyTuple_New(3);
         // ... get operator size
         pValue = PyLong_FromLong(this->mspSetup->specSize());
         PyTuple_SetItem(pArgs, 0, pValue);
         // ... create boundray condition (none)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
         PyTuple_SetItem(pArgs, 1, pValue);
         // ... set coefficient to scale factor
         pValue = PyFloat_FromDouble(this->mCScale);
         PyTuple_SetItem(pArgs, 2, pValue);

         // Call d1
         PythonWrapper::setFunction("d1");
         pValue = PythonWrapper::callFunction(pArgs);
      #elif defined GEOMHDISCC_TRANSOP_BACKWARD
         // Prepare arguments to i1(...) call
         PyObject *pArgs, *pValue;
         pArgs = PyTuple_New(3);
         // ... get operator size
         pValue = PyLong_FromLong(this->mspSetup->specSize());
         PyTuple_SetItem(pArgs, 0, pValue);
         // ... create boundray condition (last mode is zero)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(99));
         PyTuple_SetItem(pArgs, 1, pValue);
         // ... set coefficient to 1.0
         pValue = PyFloat_FromDouble(1.0);
         PyTuple_SetItem(pArgs, 2, pValue);

         // Call i1
         PythonWrapper::setFunction("i1");
         pValue = PythonWrapper::callFunction(pArgs);
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Fill matrix and cleanup
      PythonWrapper::fillMatrix(this->mDiff, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
         // Factorize matrix and free memory
         this->mSDiff.compute(this->mDiff);
         // Check for successful factorisation
         if(this->mSDiff.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward differentiation failed!");
         }
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD
   }
   #endif //defined GEOMHDISCC_TRANSOP_FORWARD || defined GEOMHDISCC_TRANSOP_BACKWARD

   #if defined GEOMHDISCC_TRANSOP_RECURRENCE
   void ChebyshevFftwTransform::recurrenceDiff(Matrix& rDealiased, const Matrix& chebVal) const
   {
      int i = chebVal.rows()-1;

      // Set T_N to zero
      rDealiased.row(i).setConstant(0.0);
      --i;

      // Compute T_N-1
      rDealiased.row(i) = static_cast<MHDFloat>(2*(i+1))*chebVal.row(i+1);
      --i;

      // Compute remaining modes
      for(; i >= 0; --i)
      {
         rDealiased.row(i) = rDealiased.row(i+2) + static_cast<MHDFloat>(2*(i+1))*chebVal.row(i+1);
      }
   }

   void ChebyshevFftwTransform::recurrenceDiff(Matrix& rDealiased, const MatrixZ& chebVal, const bool useImag) const
   {
      int i = chebVal.rows()-1;

      // Set T_N to zero
      rDealiased.row(i).setConstant(0.0);
      --i;

      if(useImag)
      {
         // Compute T_N-1
         rDealiased.row(i) = static_cast<MHDFloat>(2*(i+1))*chebVal.row(i+1).imag();
          --i;

         // Compute remaining modes
         for(; i >= 0; --i)
         {
            rDealiased.row(i) = rDealiased.row(i+2) + static_cast<MHDFloat>(2*(i+1))*chebVal.row(i+1).imag();
         }
      } else
      {
         // Compute T_N-1
         rDealiased.row(i) = static_cast<MHDFloat>(2*(i+1))*chebVal.row(i+1).real();
         --i;

         // Compute remaining modes
         for(; i >= 0; --i)
         {
            rDealiased.row(i) = rDealiased.row(i+2) + static_cast<MHDFloat>(2*(i+1))*chebVal.row(i+1).real();
         }
      }
   }
   #endif //defined GEOMHDISCC_TRANSOP_RECURRENCE

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
