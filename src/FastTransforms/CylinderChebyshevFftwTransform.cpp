/** 
 * @file CylinderChebyshevFftwTransform.cpp
 * @brief Source of the implementation of the Chebyshev FFTW transform for a cylinder radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/CylinderChebyshevFftwTransform.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array CylinderChebyshevFftwTransform::generateGrid(const int size)
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

   CylinderChebyshevFftwTransform::CylinderChebyshevFftwTransform()
      : mFPlan(NULL), mBPlan(NULL)
   {
   }

   CylinderChebyshevFftwTransform::~CylinderChebyshevFftwTransform()
   {
      // Cleanup memory used by FFTW
      FftwLibrary::cleanupFft();
   }

   void CylinderChebyshevFftwTransform::init(CylinderChebyshevFftwTransform::SharedSetupType spSetup)
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

   void CylinderChebyshevFftwTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void CylinderChebyshevFftwTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      //
      // No possible options
      //
   }

   Array CylinderChebyshevFftwTransform::meshGrid() const
   {
      return CylinderChebyshevFftwTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void CylinderChebyshevFftwTransform::initFft()
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

   void CylinderChebyshevFftwTransform::initOperators()
   {
      // First derivative
      this->mDiffE.resize(this->mspSetup->specSize(),this->mspSetup->specSize());
      this->mDiffO.resize(this->mspSetup->specSize(),this->mspSetup->specSize());
      // Division by R
      this->mDivRE.resize(this->mspSetup->specSize(),this->mspSetup->specSize());
      this->mDivRO.resize(this->mspSetup->specSize(),this->mspSetup->specSize());

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.cylindrical.cylinder_radius");

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
         PyTuple_SetItem(pArgs, 2, pValue);

         // ... set even parity
         pValue = PyLong_FromLong(0);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call d1
         PythonWrapper::setFunction("d1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix and clenup
         PythonWrapper::fillMatrix(this->mDiffE, pValue);
         Py_DECREF(pValue);

         // ... set odd parity
         pValue = PyLong_FromLong(1);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call d1
         PythonWrapper::setFunction("d1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix
         PythonWrapper::fillMatrix(this->mDiffO, pValue);
         Py_DECREF(pValue);
      #elif defined GEOMHDISCC_TRANSOP_BACKWARD
         // Prepare arguments to i1(...) call
         PyObject *pArgs, *pValue;
         pArgs = PyTuple_New(3);
         // ... get operator size
         pValue = PyLong_FromLong(this->mspSetup->specSize());
         PyTuple_SetItem(pArgs, 0, pValue);
         // ... create boundray condition (last mode zero)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(99));
         PyTuple_SetItem(pArgs, 2, pValue);

         // ... set even parity
         pValue = PyLong_FromLong(0);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call i1
         PythonWrapper::setFunction("i1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix and clenup
         PythonWrapper::fillMatrix(this->mDiffE, pValue);
         Py_DECREF(pValue);

         // ... set odd parity
         pValue = PyLong_FromLong(1);
         PyTuple_SetItem(pArgs, 1, pValue);

         // Call i1
         PythonWrapper::setFunction("i1");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix
         PythonWrapper::fillMatrix(this->mDiffO, pValue);
         Py_DECREF(pValue);
      #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Fill matrix and cleanup
      PythonWrapper::finalize();

      #if defined GEOMHDISCC_TRANSOP_BACKWARD
         // Factorize matrix and free memory
         this->mSDiffE.compute(this->mDiffE);
         // Check for successful factorisation
         if(this->mSDiffE.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward even differentiation failed!");
         }
         this->mDiffE.resize(0,0);
         // Factorize matrix and free memory
         this->mSDiffO.compute(this->mDiffO);
         // Check for successful factorisation
         if(this->mSDiffO.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward odd differentiation failed!");
         }
         this->mDiffO.resize(0,0);

         // Factorize division matrix and free memory
         this->mSDivRE.compute(this->mDivRE);
         // Check for successful factorisation
         if(this->mSDivRE.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward even division failed!");
         }
         this->mDivRE.resize(0,0);
         // Factorize division matrix and free memory
         this->mSDivRO.compute(this->mDivRO);
         // Check for successful factorisation
         if(this->mSDivRO.info() != Eigen::Success)
         {
            throw Exception("Factorization of backward odd division failed!");
         }
         this->mDivRO.resize(0,0);
      #endif //defined GEOMHDISCC_TRANSOP_BACKWARD
   }

   void CylinderChebyshevFftwTransform::cleanupFft()
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
   MHDFloat CylinderChebyshevFftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
