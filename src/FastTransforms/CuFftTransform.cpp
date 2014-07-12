/** 
 * @file CuFftTransform.cpp
 * @brief Source of the implementation of the cuFFT transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/CuFftTransform.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "FastTransforms/CuFftLibrary.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array CuFftTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create equispaced FFT grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = 2.0*Math::PI*static_cast<MHDFloat>(k)/static_cast<MHDFloat>(size);
      }

      return grid;
   }

   CuFftTransform::CuFftTransform()
      : mFPlan(), mBPlan()
   {
   }

   CuFftTransform::~CuFftTransform()
   {  
      // Cleanup memory used by cuFFT
      CuFftLibrary::cleanupFft();
   }

   void CuFftTransform::init(CuFftTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Set the scaling factor
      this->mspSetup->setScale(1.0/static_cast<MHDFloat>(this->mspSetup->fwdSize()));

      // Initialise cuFFT plans
      this->initFft();

      // Register the cuFFT object
      CuFftLibrary::registerFft();
   }

   void CuFftTransform::requiredOptions(std::set<NonDimensional::Id>& list) const
   {
      //
      // No possible options
      //
   }

   void CuFftTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options)
   {
      //
      // No possible options
      //
   }

   Array CuFftTransform::meshGrid() const
   {
      return CuFftTransform::generateGrid(this->mspSetup->fwdSize());
   }

   void CuFftTransform::initFft()
   {
      int fwdSize = this->mspSetup->fwdSize();
      int bwdSize = this->mspSetup->bwdSize();
      int howmany = this->mspSetup->howmany();

      // Create the two plans
      int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::MIXED)
      {
         // Create the real to complex plan
         checkCudaErrors(cufftPlanMany(&this->mFPlan, 1, fftSize, NULL, 1, fwdSize, NULL, 1, bwdSize, CUFFT_D2Z, howmany));

         // Create the complex to real plan
         checkCudaErrors(cufftPlanMany(&this->mBPlan,1, fftSize, NULL, 1, bwdSize, NULL, 1, fwdSize, CUFFT_Z2D, howmany));

         // Allocate common device memory
         checkCudaErrors(cudaMalloc((void **)&this->mDevR, sizeof(cufftDoubleReal)*fwdSize*howmany));
         checkCudaErrors(cudaMalloc((void **)&this->mDevZI, sizeof(cufftDoubleComplex)*bwdSize*howmany));

      } else
      {
         // Create the forward complex to complex plan
         checkCudaErrors(cufftPlanMany(&this->mFPlan, 1, fftSize, NULL, 1, fwdSize, NULL, 1, bwdSize, CUFFT_Z2Z, howmany));

         // Create the backward complex to complex plan
         checkCudaErrors(cufftPlanMany(&this->mBPlan, 1, fftSize, NULL, 1, bwdSize, NULL, 1, fwdSize, CUFFT_Z2Z, howmany));

         // Allocate common device memory
         checkCudaErrors(cudaMalloc((void **)&this->mDevZI, sizeof(cufftDoubleComplex)*bwdSize*howmany));
         checkCudaErrors(cudaMalloc((void **)&this->mDevZO, sizeof(cufftDoubleComplex)*fwdSize*howmany));
      }

      // Initialise temporary storage
      this->mTmpRIn.setZero(bwdSize, howmany);
      this->mTmpZIn.setZero(fwdSize, howmany);
   }

   void CuFftTransform::cleanupFft()
   {
      // Detroy forward plan
      if(this->mFPlan)
      {
         checkCudaErrors(cufftDestroy(this->mFPlan));
      }

      // Detroy backward plan
      if(this->mBPlan)
      {
         checkCudaErrors(cufftDestroy(this->mBPlan));
      }

      // Unregister the cuFFT object
      CuFftLibrary::unregisterFft();

      // cleanup cuFFT library
      CuFftLibrary::cleanupFft();
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat CuFftTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the cuFFT plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
