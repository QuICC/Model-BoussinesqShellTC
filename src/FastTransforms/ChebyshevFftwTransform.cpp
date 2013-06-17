/** \file ChebyshevFftwTransform.cpp
 *  \brief Source of the implementation of the FFTW transform
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

namespace GeoMHDiSCC {

namespace Transform {

   Array ChebyshevFftwTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create Chebyshev grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = std::cos((MathConstants::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(size));
      }

      return grid;
   }

   ChebyshevFftwTransform::ChebyshevFftwTransform()
      : mFPlan(NULL), mBPlan(NULL)
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

      // Register the FFTW object
      FftwLibrary::registerFft();
   }

   Array ChebyshevFftwTransform::meshGrid() const
   {
      return ChebyshevFftwTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void ChebyshevFftwTransform::initFft()
   {  
      /// \mhdBug implement strideed stranforms for complex <-> complex case if possible

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
