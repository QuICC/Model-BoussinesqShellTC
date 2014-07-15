/**
 * @file CuFftTransform.hpp
 * @brief Implementation of the FFTW transform 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CUFFTTRANSFORM_HPP
#define CUFFTTRANSFORM_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//

// System includes
//
#include <set>
#include <map>

// External includes
//
#include <cuda_runtime.h>
#include <cufft.h>

#ifdef GEOMHDISCC_DEBUG
   #include <helper_cuda.h>
#else
   #define checkCudaErrors(x) x
#endif //GEOMHDISCC_DEBUG

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "Enums/Arithmetics.hpp"
#include "FastTransforms/FftSetup.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Simple struct holding details about FFT transform
    */
   struct FftIds {

      /**
       * @brief Simple struct holding the projector IDs
       */
      struct Projectors
      {
         /// Enum of projector IDs
         enum Id {PROJ,DIFF};
      };

      /**
       * @brief Simple struct holding the integrator IDs
       */
      struct Integrators
      {
         /// Enum of integrator IDs
         enum Id {INTG};
      };

   };

   /**
    * @brief Interface class to the FFTW routines
    */ 
   class CuFftTransform 
   {
      public:
         /// Typedef for the configuration class
         typedef FftSetup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedFftSetup SharedSetupType;

         /// Typedef for the Projector type
         typedef FftIds::Projectors ProjectorType;

         /// Typedef for the Integrator type
         typedef FftIds::Integrators IntegratorType;

         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size); 

         /**
          * @brief Very basic constructor
          */
         CuFftTransform();

         /**
          * @brief Destroy the FFTW plans
          */
         ~CuFftTransform();

         /**
          * @brief Initialise the FFT computations (plans, etc)
          *
          * Compute the optimal plans for the required FFT transforms
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief set list of required options
          */
         void requiredOptions(std::set<NonDimensional::Id>& list) const;

         /**
          * @brief Set the required options
          */
         void setOptions(const std::map<NonDimensional::Id, MHDFloat>& options);

         /**
          * @brief Get the physical grid
          */
         Array meshGrid() const; 

         /**
          * @brief Compute forward FFT (R2C)
          *
          * Compute the FFT from real physical space to complex spectral space
          *
          * @param rFFTVal Output FFT transformed values
          * @param physVal Input physical values
          * @param integrator Integrator to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void integrate(MatrixZ& rFFTVal, const Matrix& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute backward FFT (C2R)
          *
          * Compute the FFT from complex spectral space to real physical space
          *
          * @param rPhysVal Output physical values
          * @param fftVal Input FFT values
          * @param projector  Projector to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void project(Matrix& rPhysVal, const MatrixZ& fftVal, ProjectorType::Id projector);

         /**
          * @brief Compute forward FFT (C2C)
          *
          * Compute the FFT from complex physical space to complex spectral space
          *
          * @param rFFTVal Output FFT transformed values
          * @param physVal Input physical values
          * @param integrator Integrator to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void integrate(MatrixZ& rFFTVal, const MatrixZ& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute backward FFT (C2C)
          *
          * Compute the FFT from complex spectral space to complex physical space
          *
          * @param rPhysVal   Output physical values
          * @param fftVal     Input FFT values
          * @param projector  Projector to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void project(MatrixZ& rPhysVal, const MatrixZ& fftVal, ProjectorType::Id projector);

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief FFT setup object providing the sizes
          */
         SharedFftSetup    mspSetup;

         /**
          * @brief Plan for the forward transform (real->complex or complex->complex)
          */
         cufftHandle   mFPlan;

         /**
          * @brief Plan for the backward transform (complex->real or complex->complex)
          */
         cufftHandle   mBPlan;

         /**
          * @brief Temporary storage used in the projections (complex -> real)
          */
         MatrixZ  mTmpRIn;

         /**
          * @brief Temporary storage used in the projections (complex -> complex)
          */
         MatrixZ  mTmpZIn;

         /**
          * @brief Real storage on device
          */
         cufftDoubleReal *mDevR;

         /**
          * @brief Complex storage on device
          */
         cufftDoubleComplex *mDevZI;

         /**
          * @brief Complex storage on device
          */
         cufftDoubleComplex *mDevZO;

         /**
          * @brief Initialise the FFTW transforms (i.e. create plans, etc)
          */
         void initFft();

         /**
          * @brief Cleanup memory used by FFTW on destruction
          */
         void cleanupFft();
   };

   template <Arithmetics::Id TOperation> void CuFftTransform::integrate(MatrixZ& rFFTVal, const Matrix& physVal, CuFftTransform::IntegratorType::Id integrator)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::MIXED);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rFFTVal.rows() == this->mspSetup->bwdSize());
      assert(rFFTVal.cols() == this->mspSetup->howmany());

      // Do transform
      checkCudaErrors(cudaMemcpy(this->mDevR, physVal.data(), sizeof(cufftDoubleReal)*physVal.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecD2Z(this->mFPlan, this->mDevR, this->mDevZI));
      checkCudaErrors(cudaMemcpy(rFFTVal.data(), this->mDevZI, sizeof(cufftDoubleComplex)*rFFTVal.size(), cudaMemcpyDeviceToHost));

      // Rescale output from FFT
      rFFTVal *= this->mspSetup->scale();
   }

   template <Arithmetics::Id TOperation> void CuFftTransform::project(Matrix& rPhysVal, const MatrixZ& fftVal, CuFftTransform::ProjectorType::Id projector)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

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
      checkCudaErrors(cudaMemcpy(this->mDevZI, this->mTmpRIn.data(), sizeof(cufftDoubleComplex)*this->mTmpRIn.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecZ2D(this->mBPlan, this->mDevZI, this->mDevR));
      checkCudaErrors(cudaMemcpy(rPhysVal.data(), this->mDevR, sizeof(cufftDoubleReal)*rPhysVal.size(), cudaMemcpyDeviceToHost));
   }

   template <Arithmetics::Id TOperation> void CuFftTransform::integrate(MatrixZ& rFFTVal, const MatrixZ& physVal, CuFftTransform::IntegratorType::Id integrator)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // Assert that a non mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPLEX);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rFFTVal.rows() == this->mspSetup->bwdSize());
      assert(rFFTVal.cols() == this->mspSetup->howmany());

      // Do transform
      checkCudaErrors(cudaMemcpy(this->mDevZI, physVal.data(), sizeof(cufftDoubleComplex)*physVal.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecZ2Z(this->mFPlan, this->mDevZI, this->mDevZO, CUFFT_FORWARD));
      checkCudaErrors(cudaMemcpy(rFFTVal.data(), this->mDevZO, sizeof(cufftDoubleComplex)*rFFTVal.size(), cudaMemcpyDeviceToHost));

      // Rescale output from FFT
      rFFTVal *= this->mspSetup->scale();
   }

   template <Arithmetics::Id TOperation> void CuFftTransform::project(MatrixZ& rPhysVal, const MatrixZ& fftVal, CuFftTransform::ProjectorType::Id projector)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

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
      checkCudaErrors(cudaMemcpy(this->mDevZI, this->mTmpZIn.data(), sizeof(cufftDoubleComplex)*this->mTmpZIn.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecZ2Z(this->mFPlan, this->mDevZI, this->mDevZO, CUFFT_INVERSE));
      checkCudaErrors(cudaMemcpy(rPhysVal.data(), this->mDevZO, sizeof(cufftDoubleComplex)*rPhysVal.size(), cudaMemcpyDeviceToHost));
   }

}
}

#endif // CUFFTTRANSFORM_HPP
