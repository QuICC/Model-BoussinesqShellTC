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

      // Register the cuFFT object
      CuFftLibrary::registerFft();

      // Initialise cuFFT plans
      this->initFft();
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

      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         this->mStreamBatch.push_back(howmany/CuFftLibrary::NSTREAMS);
         if(i == CuFftLibrary::NSTREAMS - 1)
         {
            this->mStreamBatch.back() += howmany % CuFftLibrary::NSTREAMS;
         }
      }

      // Create the two plans
      int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::MIXED)
      {
         for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
         {
            // Initialise handles
            this->mFPlan.push_back(cufftHandle());
            this->mBPlan.push_back(cufftHandle());

            // Create the real to complex plan
            checkCudaErrors(cufftPlanMany(&this->mFPlan.at(i), 1, fftSize, NULL, 1, fwdSize, NULL, 1, bwdSize, CUFFT_D2Z, this->mStreamBatch.at(i)));
            checkCudaErrors(cufftSetStream(this->mFPlan.at(i), CuFftLibrary::sStream.at(i)));

            // Create the complex to real plan
            checkCudaErrors(cufftPlanMany(&this->mBPlan.at(i),1, fftSize, NULL, 1, bwdSize, NULL, 1, fwdSize, CUFFT_Z2D, this->mStreamBatch.at(i)));
            checkCudaErrors(cufftSetStream(this->mBPlan.at(i), CuFftLibrary::sStream.at(i)));

            // Allocate common device memory
            cufftDoubleReal * tmpR;
            checkCudaErrors(cudaMalloc((void **)&tmpR, sizeof(cufftDoubleReal)*fwdSize*this->mStreamBatch.at(i)));
            this->mpDevR.push_back(tmpR);
            cufftDoubleComplex * tmpZ;
            checkCudaErrors(cudaMalloc((void **)&tmpZ, sizeof(cufftDoubleComplex)*bwdSize*this->mStreamBatch.at(i)));
            this->mpDevZI.push_back(tmpZ);
         }

      } else
      {
         for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
         {
            // Initialise handles
            this->mFPlan.push_back(cufftHandle());
            this->mBPlan.push_back(cufftHandle());

            // Create the forward complex to complex plan
            checkCudaErrors(cufftPlanMany(&this->mFPlan.back(), 1, fftSize, NULL, 1, fwdSize, NULL, 1, bwdSize, CUFFT_Z2Z, this->mStreamBatch.at(i)));

            // Create the backward complex to complex plan
            checkCudaErrors(cufftPlanMany(&this->mBPlan.back(), 1, fftSize, NULL, 1, bwdSize, NULL, 1, fwdSize, CUFFT_Z2Z, this->mStreamBatch.at(i)));

            // Allocate common device memory
            cufftDoubleComplex * tmpZ;
            checkCudaErrors(cudaMalloc((void **)&tmpZ, sizeof(cufftDoubleComplex)*bwdSize*this->mStreamBatch.at(i)));
            this->mpDevZI.push_back(tmpZ);
            checkCudaErrors(cudaMalloc((void **)&tmpZ, sizeof(cufftDoubleComplex)*fwdSize*this->mStreamBatch.at(i)));
            this->mpDevZO.push_back(tmpZ);
         }
      }

      // Initialise temporary storage
      this->mTmpRIn.setZero(bwdSize, howmany);
      this->mTmpZIn.setZero(fwdSize, howmany);
   }

   void CuFftTransform::cleanupFft()
   {
      for(size_t i = 0; i < this->mFPlan.size(); i++)
      {
         // Destroy forward plan
         if(this->mFPlan.at(i))
         {
            checkCudaErrors(cufftDestroy(this->mFPlan.at(i)));
         }
      }

      for(size_t i = 0; i < this->mBPlan.size(); i++)
      {
         // Destroy backward plan
         if(this->mBPlan.at(i))
         {
            checkCudaErrors(cufftDestroy(this->mBPlan.at(i)));
         }
      }

      for(size_t i = 0; i < this->mpDevR.size(); i++)
      {
         cudaFree(this->mpDevR.at(i));
      }

      for(size_t i = 0; i < this->mpDevZI.size(); i++)
      {
         cudaFree(this->mpDevZI.at(i));
      }

      for(size_t i = 0; i < this->mpDevZO.size(); i++)
      {
         cudaFree(this->mpDevZO.at(i));
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
