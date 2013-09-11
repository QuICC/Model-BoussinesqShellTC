/** 
 * @file FftwTransform.cpp
 * @brief Source of the implementation of the FFTW transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/FftwTransform.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array FftwTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create equispaced FFT grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = 2.0*MathConstants::PI*static_cast<MHDFloat>(k)/static_cast<MHDFloat>(size);
      }

      return grid;
   }

   FftwTransform::FftwTransform()
      : mFPlan(NULL), mBPlan(NULL)
   {
   }

   FftwTransform::~FftwTransform()
   {
      // Cleanup memory used by FFTW
      FftwLibrary::cleanupFft();
   }

   void FftwTransform::init(FftwTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Set the scaling factor
      this->mspSetup->setScale(1.0/static_cast<MHDFloat>(this->mspSetup->fwdSize()));

      // Initialise FFTW plans
      this->initFft();

      // Register the FFTW object
      FftwLibrary::registerFft();
   }

   void FftwTransform::requiredOptions(std::set<NonDimensional::Id>& list) const
   {
      //
      // No possible options
      //
   }

   void FftwTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options)
   {
      //
      // No possible options
      //
   }

   Array FftwTransform::meshGrid() const
   {
      return FftwTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void FftwTransform::initFft()
   {
      int fwdSize = this->mspSetup->fwdSize();
      int bwdSize = this->mspSetup->bwdSize();
      int howmany = this->mspSetup->howmany();

      // Create the two plans
      const int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::MIXED)
      {
         // create temporary storage for plan computation
         Matrix    tmpReal(fwdSize, howmany);
         MatrixZ   tmpCplx(bwdSize, howmany);

         // Create the real to complex plan
         this->mFPlan = fftw_plan_many_dft_r2c(1, fftSize, howmany, tmpReal.data(), NULL, 1, fwdSize, reinterpret_cast<fftw_complex* >(tmpCplx.data()), NULL, 1, bwdSize, FftwLibrary::planFlag());

         // Create the complex to real plan
         this->mBPlan = fftw_plan_many_dft_c2r(1, fftSize, howmany, reinterpret_cast<fftw_complex* >(tmpCplx.data()), NULL, 1, bwdSize, tmpReal.data(), NULL, 1, fwdSize, FftwLibrary::planFlag());
      } else
      {
         MatrixZ   tmpCplxA(fwdSize, howmany);
         MatrixZ   tmpCplxB(bwdSize, howmany);

         // Create the forward complex to complex plan
         this->mFPlan = fftw_plan_many_dft(1, fftSize, howmany, reinterpret_cast<fftw_complex* >(tmpCplxA.data()), NULL, 1, fwdSize, reinterpret_cast<fftw_complex* >(tmpCplxB.data()), NULL, 1, bwdSize, FFTW_FORWARD, FftwLibrary::planFlag());

         // Create the backward complex to complex plan
         this->mBPlan = fftw_plan_many_dft(1, fftSize, howmany, reinterpret_cast<fftw_complex* >(tmpCplxB.data()), NULL, 1, bwdSize, reinterpret_cast<fftw_complex* >(tmpCplxA.data()), NULL, 1, fwdSize, FFTW_BACKWARD, FftwLibrary::planFlag());
      }

      // Initialise temporary storage
      this->mTmpRIn.setZero(bwdSize, howmany);
      this->mTmpZIn.setZero(fwdSize, howmany);
   }

   void FftwTransform::cleanupFft()
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

      // cleanup FFTW library
      FftwLibrary::cleanupFft();
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat FftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
