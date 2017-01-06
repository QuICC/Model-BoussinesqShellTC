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

   void CuFftTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void CuFftTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
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

   void CuFftTransform::integrate(MatrixZ& rFFTVal, const Matrix& physVal, CuFftTransform::IntegratorType::Id integrator)
   {
      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::MIXED);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rFFTVal.rows() == this->mspSetup->bwdSize());
      assert(rFFTVal.cols() == this->mspSetup->howmany());

      // Do transform
      int offset = 0;
      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         int sze = physVal.rows()*this->mStreamBatch.at(i);

         checkCudaErrors(cudaMemcpyAsync(this->mpDevR.at(i), physVal.data() + offset, sizeof(cufftDoubleReal)*sze, cudaMemcpyHostToDevice, CuFftLibrary::sStream.at(i)));

         offset += sze;
      }

      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         checkCudaErrors(cufftExecD2Z(this->mFPlan.at(i), this->mpDevR.at(i), this->mpDevZI.at(i)));
      }

      offset = 0;
      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         int sze = rFFTVal.rows()*this->mStreamBatch.at(i);

         checkCudaErrors(cudaMemcpyAsync(rFFTVal.data() + offset, this->mpDevZI.at(i), sizeof(cufftDoubleComplex)*sze, cudaMemcpyDeviceToHost, CuFftLibrary::sStream.at(i)));

         offset += sze;
      }

      // Rescale output from FFT
      rFFTVal *= this->mspSetup->scale();
   }

   void CuFftTransform::project(Matrix& rPhysVal, const MatrixZ& fftVal, CuFftTransform::ProjectorType::Id projector)
   {
      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::MIXED);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input  matrix
      assert(fftVal.rows() == this->mspSetup->bwdSize());
      assert(fftVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative
      if(projector == CuFftTransform::ProjectorType::DIFF)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);

         // Rescale results
         this->mTmpRIn.topRows(this->mspSetup->specSize()) = factor.asDiagonal()*fftVal.topRows(this->mspSetup->specSize());

      // Compute simple projection
      } else
      {
         // Rescale results
         this->mTmpRIn.topRows(this->mspSetup->specSize()) = fftVal.topRows(this->mspSetup->specSize());
      }

      // Set the m=0 values to zero
      this->mTmpRIn.row(0).imag().setConstant(0);

      // Set the padded values to zero
      this->mTmpRIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform
      int offset = 0;
      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         int sze = this->mTmpRIn.rows()*this->mStreamBatch.at(i);

         checkCudaErrors(cudaMemcpyAsync(this->mpDevZI.at(i), this->mTmpRIn.data() + offset, sizeof(cufftDoubleComplex)*sze, cudaMemcpyHostToDevice, CuFftLibrary::sStream.at(i)));

         offset += sze;
      }

      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         checkCudaErrors(cufftExecZ2D(this->mBPlan.at(i), this->mpDevZI.at(i), this->mpDevR.at(i)));
      }

      offset = 0;
      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         int sze = rPhysVal.rows()*this->mStreamBatch.at(i);

         checkCudaErrors(cudaMemcpyAsync(rPhysVal.data() + offset, this->mpDevR.at(i), sizeof(cufftDoubleReal)*sze, cudaMemcpyDeviceToHost, CuFftLibrary::sStream.at(i)));

         offset += sze;
      }
   }

   void CuFftTransform::integrate(MatrixZ& rFFTVal, const MatrixZ& physVal, CuFftTransform::IntegratorType::Id integrator)
   {
      // Assert that a non mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPLEX);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rFFTVal.rows() == this->mspSetup->bwdSize());
      assert(rFFTVal.cols() == this->mspSetup->howmany());

      // Do transform
      int offset = 0;
      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         int sze = physVal.rows()*this->mStreamBatch.at(i);

         checkCudaErrors(cudaMemcpyAsync(this->mpDevZI.at(i), physVal.data() + offset, sizeof(cufftDoubleComplex)*sze, cudaMemcpyHostToDevice, CuFftLibrary::sStream.at(i)));

         offset += sze;
      }
      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         checkCudaErrors(cufftExecZ2Z(this->mFPlan.at(i), this->mpDevZI.at(i), this->mpDevZO.at(i), CUFFT_FORWARD));
      }

      offset = 0;
      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         int sze = rFFTVal.rows()*this->mStreamBatch.at(i);

         checkCudaErrors(cudaMemcpyAsync(rFFTVal.data() + offset, this->mpDevZO.at(i), sizeof(cufftDoubleComplex)*sze, cudaMemcpyDeviceToHost, CuFftLibrary::sStream.at(i)));

         offset += sze;
      }

      // Rescale output from FFT
      rFFTVal *= this->mspSetup->scale();
   }

   void CuFftTransform::project(MatrixZ& rPhysVal, const MatrixZ& fftVal, CuFftTransform::ProjectorType::Id projector)
   {
      // Assert that a non mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPLEX);

      // assert on the padding size
      assert(this->mspSetup->padSize() >= 0);
      assert(this->mspSetup->bwdSize() - this->mspSetup->padSize() >= 0);

      // assert right sizes for input  matrix
      assert(fftVal.rows() == this->mspSetup->bwdSize());
      assert(fftVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Get size of positive and negative frequency parts
      int negN = this->mspSetup->specSize()/2;
      int posN = negN + (this->mspSetup->specSize()%2);

      // Compute first derivative
      if(projector == CuFftTransform::ProjectorType::DIFF)
      {
         // Get differentiation factors
         ArrayZ factor = this->mspSetup->boxScale()*Math::cI*Array::LinSpaced(posN, 0, posN-1);
         ArrayZ rfactor = this->mspSetup->boxScale()*Math::cI*(Array::LinSpaced(negN, 0, negN-1).array() - static_cast<MHDFloat>(negN));

         // Split positive and negative frequencies and compute derivative
         this->mTmpZIn.topRows(posN) = factor.asDiagonal()*fftVal.topRows(posN);
         this->mTmpZIn.bottomRows(negN) = rfactor.asDiagonal()*fftVal.bottomRows(negN);

      // Compute simple projection
      } else
      {
         // Split positive and negative frequencies
         this->mTmpZIn.topRows(posN) = fftVal.topRows(posN);
         this->mTmpZIn.bottomRows(negN) = fftVal.bottomRows(negN);
      }

      // Set the padded values to zero
      this->mTmpZIn.block(posN, 0, this->mspSetup->padSize(), this->mTmpZIn.cols()).setZero();

      // Do transform
      int offset = 0;
      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         int sze = this->mTmpZIn.rows()*this->mStreamBatch.at(i);

         checkCudaErrors(cudaMemcpyAsync(this->mpDevZI.at(i), this->mTmpZIn.data()+offset, sizeof(cufftDoubleComplex)*sze, cudaMemcpyHostToDevice, CuFftLibrary::sStream.at(i)));

         offset += sze;
      }

      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         checkCudaErrors(cufftExecZ2Z(this->mFPlan.at(i), this->mpDevZI.at(i), this->mpDevZO.at(i), CUFFT_INVERSE));
      }

      offset = 0;
      for(int i = 0; i < CuFftLibrary::NSTREAMS; i++)
      {
         int sze = rPhysVal.rows()*this->mStreamBatch.at(i);

         checkCudaErrors(cudaMemcpyAsync(rPhysVal.data()+offset, this->mpDevZO.at(i), sizeof(cufftDoubleComplex)*sze, cudaMemcpyDeviceToHost, CuFftLibrary::sStream.at(i)));

         offset += sze;
      }
   }

#ifdef QUICC_STORAGEPROFILE
   MHDFloat CuFftTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the cuFFT plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}
