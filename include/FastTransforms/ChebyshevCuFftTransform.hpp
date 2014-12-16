/**
 * @file ChebyshevCuFftTransform.hpp
 * @brief Implementation of the cuFFT transform for a Chebyshev expansion 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CHEBYSHEVCUFFTTRANSFORM_HPP
#define CHEBYSHEVCUFFTTRANSFORM_HPP

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
#include "Enums/Dimensions.hpp"
#include "Enums/Arithmetics.hpp"
#include "Enums/NonDimensional.hpp"
#include "FastTransforms/FftSetup.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "SparseSolvers/SparseLinearSolverTools.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Simple struct holding details about ChebyshevFFT transform
    */
   struct ChebyshevFftIds {

      /**
       * @brief Simple struct holding the projector IDs
       */
      struct Projectors
      {
         /** Enum of projector IDs
          *    - PROJ: projection
          *    - DIFF: D
          *    - DIVR: 1/r
          *    - DIVR2: 1/r^2
          *    - DIVRDIFFR: 1/r D r
          *    - DIFFDIVR: D 1/r
          */
         enum Id {PROJ,  DIFF, DIVR, DIVR2, DIVRDIFFR, DIFFDIVR};
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
    * @brief Implementation of the cuFFT transform for a Chebyshev expansion
    */ 
   class ChebyshevCuFftTransform
   {
      public:
         /// Typedef for the configuration class
         typedef FftSetup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedFftSetup SharedSetupType;

         /// Typedef for the Projector type
         typedef ChebyshevFftIds::Projectors ProjectorType;

         /// Typedef for the Integrator type
         typedef ChebyshevFftIds::Integrators IntegratorType;

         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size); 

         /**
          * @brief Very basic constructor
          */
         ChebyshevCuFftTransform();

         /**
          * @brief Destroy the cuFFT plans
          */
         ~ChebyshevCuFftTransform();

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
         void requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const;

         /**
          * @brief Set the required options
          */
         void setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId);

         /**
          * @brief Get the physical grid
          */
         Array meshGrid() const; 

         /**
          * @brief Compute forward FFT (R2R)
          *
          * Compute the FFT from real physical space to Chebyshev spectral space
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void integrate(Matrix& rChebVal, const Matrix& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute forward FFT (C2C)
          *
          * Compute the FFT from real physical space to Chebyshev spectral space
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void integrate(MatrixZ& rChebVal, const MatrixZ& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute backward FFT (R2R)
          *
          * Compute the FFT from Chebyshev spectral space to real physical space
          *
          * @param rPhysVal   Output physical values
          * @param chebVal    Input Chebyshev coefficients
          * @param projector  Projector to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void project(Matrix& rPhysVal, const Matrix& chebVal, ProjectorType::Id projector);

         /**
          * @brief Compute backward FFT (C2C)
          *
          * Compute the FFT from Chebyshev spectral space to real physical space
          *
          * @param rPhysVal   Output physical values
          * @param chebVal    Input Chebyshev coefficients
          * @param projector  Projector to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ProjectorType::Id projector);

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
          * @brief Plan for the forward transform (real -> real)
          */
         cufftHandle   mFPlan;

         /**
          * @brief Plan for the backward transform (real -> real)
          */
         cufftHandle   mBPlan;

         /**
          * @brief Scale factor for Chebyshev axis
          */
         MHDFloat mCScale;

         /**
          * @brief Temporary real storage
          */
         Matrix   mTmpR;

         /**
          * @brief Temporary complex storage
          */
         MatrixZ  mTmpZ;

         /**
          * @brief Real storage on device
          */
         cufftDoubleReal *mDevR;

         /**
          * @brief Complex storage on device
          */
         cufftDoubleComplex *mDevZ;

         /**
          * @brief Storage for the Chebyshev differentiation matrix
          */
         SparseMatrix   mDiff;

         /**
          * @brief Phase shift to convert DFT to DCT
          */
         ArrayZ   mPhase;

         /**
          * @brief 1/Phase shift to convert DFT to DCT
          */
         ArrayZ   mPhase_1;

         /**
          * @brief Initialise the spectral operators
          */
         void initOperators();

         /**
          * @brief Initialise the cuFFT transforms (i.e. create plans, etc)
          */
         void initFft();

         /**
          * @brief Cleanup memory used by cuFFT on destruction
          */
         void cleanupFft();
   };

   template <Arithmetics::Id TOperation> void ChebyshevCuFftTransform::integrate(Matrix& rChebVal, const Matrix& physVal, ChebyshevCuFftTransform::IntegratorType::Id integrator)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // Assert that a mixed transform was not setup
      assert(this->mspSetup->type() == FftSetup::REAL);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Copy data to perform 2N R2C
      this->mTmpR.topRows(this->mspSetup->fwdSize()) = physVal;
      this->mTmpR.bottomRows(this->mspSetup->fwdSize()) = physVal.colwise().reverse();

      // Do transform
      checkCudaErrors(cudaMemcpy(this->mDevR, this->mTmpR.data(), sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecD2Z(this->mFPlan, this->mDevR, this->mDevZ));
      checkCudaErrors(cudaMemcpy(this->mTmpZ.data(), this->mDevZ, sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyDeviceToHost));

      // Extract real part and rescale to remove FFT scaling
      this->mTmpZ = this->mPhase_1.asDiagonal()*this->mTmpZ;
      rChebVal = this->mspSetup->scale()*this->mTmpZ.topRows(this->mspSetup->bwdSize()).real();
   }

   template <Arithmetics::Id TOperation> void ChebyshevCuFftTransform::project(Matrix& rPhysVal, const Matrix& chebVal, ChebyshevCuFftTransform::ProjectorType::Id projector)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

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
      if(projector == ChebyshevCuFftTransform::ProjectorType::DIFF)
      {
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = this->mDiff*chebVal.topRows(this->mspSetup->specSize());

      // Compute simple projection
      } else
      {
         // Copy into other array
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = chebVal.topRows(this->mspSetup->specSize());
      }

      // Set the padded values to zero
      this->mTmpZ.imag().setZero();
      this->mTmpZ.bottomRows(this->mspSetup->padSize()+1).real().setZero();
      this->mTmpZ = this->mPhase.asDiagonal()*this->mTmpZ;

      // Do transform
      checkCudaErrors(cudaMemcpy(this->mDevZ, this->mTmpZ.data(), sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecZ2D(this->mBPlan, this->mDevZ, this->mDevR));
      checkCudaErrors(cudaMemcpy(this->mTmpR.data(), this->mDevR, sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyDeviceToHost));

      // Extract physical values
      rPhysVal = this->mTmpR.topRows(this->mspSetup->fwdSize());
   }

   template <Arithmetics::Id TOperation> void ChebyshevCuFftTransform::integrate(MatrixZ& rChebVal, const MatrixZ& physVal, ChebyshevCuFftTransform::IntegratorType::Id integrator)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Do transform of real part
      this->mTmpR.topRows(this->mspSetup->fwdSize()) = physVal.real().colwise().reverse();
      this->mTmpR.bottomRows(this->mspSetup->fwdSize()) = physVal.real();
      checkCudaErrors(cudaMemcpy(this->mDevR, this->mTmpR.data(), sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecD2Z(this->mFPlan, this->mDevR, this->mDevZ));
      checkCudaErrors(cudaMemcpy(this->mTmpZ.data(), this->mDevZ, sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyDeviceToHost));

      // Rescale FFT output
      this->mTmpZ = this->mPhase_1.asDiagonal()*this->mTmpZ;
      rChebVal.real() = this->mspSetup->scale()*this->mTmpZ.topRows(this->mspSetup->bwdSize()).real();

      // Do transform of imaginary part
      this->mTmpR.topRows(this->mspSetup->fwdSize()) = physVal.imag().colwise().reverse();
      this->mTmpR.bottomRows(this->mspSetup->fwdSize()) = physVal.imag();
      checkCudaErrors(cudaMemcpy(this->mDevR, this->mTmpR.data(), sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecD2Z(this->mFPlan, this->mDevR, this->mDevZ));
      checkCudaErrors(cudaMemcpy(this->mTmpZ.data(), this->mDevZ, sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyDeviceToHost));

      // Rescale FFT output
      this->mTmpZ = this->mPhase_1.asDiagonal()*this->mTmpZ;
      rChebVal.imag() = this->mspSetup->scale()*this->mTmpZ.topRows(this->mspSetup->bwdSize()).real();
   }

   template <Arithmetics::Id TOperation> void ChebyshevCuFftTransform::project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ChebyshevCuFftTransform::ProjectorType::Id projector)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

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
      if(projector == ChebyshevCuFftTransform::ProjectorType::DIFF)
      {
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).real();

      // Compute simple projection of real part
      } else
      {
         // Copy values into simple matrix
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = chebVal.topRows(this->mspSetup->specSize()).real();
      }

      // Set the padded values to zero
      this->mTmpZ.imag().setZero();
      this->mTmpZ.bottomRows(this->mspSetup->padSize()+1).real().setZero();
      this->mTmpZ = this->mPhase.asDiagonal()*this->mTmpZ;

      // Do transform of real part
      checkCudaErrors(cudaMemcpy(this->mDevZ, this->mTmpZ.data(), sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecZ2D(this->mBPlan, this->mDevZ, this->mDevR));
      checkCudaErrors(cudaMemcpy(this->mTmpR.data(), this->mDevR, sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyDeviceToHost));
      rPhysVal.real() = this->mTmpR.bottomRows(this->mspSetup->fwdSize());

      // Compute first derivative of imaginary part
      if(projector == ChebyshevCuFftTransform::ProjectorType::DIFF)
      {
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).imag();

      // Compute simple projection of imaginary part
      } else
      {
         // Rescale results
         this->mTmpZ.topRows(this->mspSetup->specSize()).real() = chebVal.topRows(this->mspSetup->specSize()).imag();
      }

      // Set the padded values to zero
      this->mTmpZ.imag().setZero();
      this->mTmpZ.bottomRows(this->mspSetup->padSize()+1).real().setZero();
      this->mTmpZ = this->mPhase.asDiagonal()*this->mTmpZ;

      // Do transform of imaginary part
      checkCudaErrors(cudaMemcpy(this->mDevZ, this->mTmpZ.data(), sizeof(cufftDoubleComplex)*this->mTmpZ.size(), cudaMemcpyHostToDevice));
      checkCudaErrors(cufftExecZ2D(this->mBPlan, this->mDevZ, this->mDevR));
      checkCudaErrors(cudaMemcpy(this->mTmpR.data(), this->mDevR, sizeof(cufftDoubleReal)*this->mTmpR.size(), cudaMemcpyDeviceToHost));
      rPhysVal.imag() = this->mTmpR.bottomRows(this->mspSetup->fwdSize());
   }

}
}

#endif // CHEBYSHEVCUFFTTRANSFORM_HPP
