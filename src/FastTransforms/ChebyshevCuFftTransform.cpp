/** 
 * @file ChebyshevCuFftTransform.cpp
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
#include "FastTransforms/ChebyshevCuFftTransform.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "FastTransforms/CuFftLibrary.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array ChebyshevCuFftTransform::generateGrid(const int size)
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

   ChebyshevCuFftTransform::ChebyshevCuFftTransform()
      : mFPlan(NULL), mBPlan(NULL)
   {
   }

   ChebyshevCuFftTransform::~ChebyshevCuFftTransform()
   {
      // Cleanup memory used by FFTW
      CuFftLibrary::cleanupFft();
   }

   void ChebyshevCuFftTransform::init(ChebyshevCuFftTransform::SharedSetupType spSetup)
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
      CuFftLibrary::registerFft();
   }

   void ChebyshevCuFftTransform::requiredOptions(std::set<NonDimensional::Id>& list) const
   {
      //
      // No possible options
      //
   }

   void ChebyshevCuFftTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options)
   {
      //
      // No possible options
      //
   }

   Array ChebyshevCuFftTransform::meshGrid() const
   {
      return ChebyshevCuFftTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void ChebyshevCuFftTransform::initFft()
   {  
      /// \mhdBug implement strided stranforms for complex <-> complex case if possible

      int fwdSize = 2*this->mspSetup->fwdSize();
      int bwdSize = this->mspSetup->bwdSize()+1;
      int howmany = this->mspSetup->howmany();

      // Create the two plans
      const int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::COMPONENT)
      {
      } else
      {
      }

      // Initialise temporary storage
      this->mTmpR.setZero(fwdSize, howmany);
      this->mTmpZ.setZero(bwdSize, howmany);

      // Create the physical to spectral plan
      this->mFPlan = fftw_plan_many_dft_r2c(1, fftSize, howmany, this->mTmpR.data(), NULL, 1, fwdSize, reinterpret_cast<fftw_complex* >(this->mTmpZ.data()), NULL, 1, bwdSize, CuFftLibrary::planFlag());

      // Create the spectral to physical plan
      this->mBPlan = fftw_plan_many_dft_c2r(1, fftSize, howmany, reinterpret_cast<fftw_complex* >(this->mTmpZ.data()), NULL, 1, bwdSize, this->mTmpR.data(), NULL, 1, fwdSize, CuFftLibrary::planFlag());


      this->mPhase = (-Math::cI*2.0*Math::PI*((this->mspSetup->fwdSize()-0.5)/(2.0*this->mspSetup->fwdSize()))*Array::LinSpaced(this->mTmpZ.rows(), 0, this->mTmpZ.rows()-1)).array().exp();
      this->mPhase_1 = 1.0/this->mPhase.array();
   }

   void ChebyshevCuFftTransform::initOperators()
   {
      this->mDiff.resize(this->mspSetup->specSize(),this->mspSetup->specSize());

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.cartesian.cartesian_1d");

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
      // ... set coefficient to 1.0
      pValue = PyFloat_FromDouble(1.0);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Call d1
      PythonWrapper::setFunction("d1");
      pValue = PythonWrapper::callFunction(pArgs);

      // Fill matrix and clenup
      PythonWrapper::fillMatrix(this->mDiff, pValue);
      PythonWrapper::finalize();
   }

   void ChebyshevCuFftTransform::cleanupFft()
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
      CuFftLibrary::unregisterFft();

      // cleanup fftw library
      CuFftLibrary::cleanupFft();
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat ChebyshevCuFftTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
