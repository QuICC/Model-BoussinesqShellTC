/** \file ChebyshevFftwTransform.hpp
 *  \brief Implementation of the FFTW transform for a Chebyshev expansion
 */

#ifndef CHEBYSHEVFFTWTRANSFORM_HPP
#define CHEBYSHEVFFTWTRANSFORM_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//

// System includes
//

// External includes
//
#include <fftw3.h>

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Arithmetics.hpp"
#include "FastTransforms/FftSetup.hpp"
#include "SpectralOperators/ChebyshevOperator.hpp"

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
         /// Enum of projector IDs
         enum Id {PROJ,  DIFF};
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
    * @brief Implementation of the FFTW transform for a Chebyshev expansion
    */ 
   class ChebyshevFftwTransform
   {
      public:
         /// Typedef for the configuration class to use
         typedef SharedFftSetup SharedSetupType;

         /// Typedef for the Projector type
         typedef ChebyshevFftIds::Projectors ProjectorType;

         /// Typedef for the Integrator type
         typedef ChebyshevFftIds::Integrators IntegratorType;

         /**
          * @brief Very basic constructor
          */
         ChebyshevFftwTransform();

         /**
          * @brief Destroy the FFTW plans
          */
         ~ChebyshevFftwTransform();

         /**
          * @brief Initialise the FFT computations (plans, etc)
          *
          * Compute the optimal plans for the required FFT transforms
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

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
          * @brief FFTW plan for the forward transform (real -> real)
          */
         fftw_plan   mFPlan;

         /**
          * @brief FFTW plan for the backward transform (real -> real)
          */
         fftw_plan   mBPlan;

         /**
          * @brief Storage for data input
          */
         Matrix   mTmpIn;

         /**
          * @brief Storage for data output
          */
         Matrix   mTmpOut;

         /**
          * @brief Initialise the FFTW transforms (i.e. create plans, etc)
          */
         void initFft();

         /**
          * @brief Cleanup memory used by FFTW on destruction
          */
         void cleanupFft();
   };

   template <Arithmetics::Id TOperation> void ChebyshevFftwTransform::integrate(Matrix& rChebVal, const Matrix& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Do transform
      fftw_execute_r2r(this->mFPlan, const_cast<MHDFloat *>(physVal.data()), rChebVal.data());

      // Rescale to remove FFT scaling
      rChebVal *= this->mspSetup->scale();
   }

   template <Arithmetics::Id TOperation> void ChebyshevFftwTransform::project(Matrix& rPhysVal, const Matrix& chebVal, ChebyshevFftwTransform::ProjectorType::Id projector)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

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
      if(projector == ChebyshevFftwTransform::ProjectorType::DIFF)
      {
         Spectral::ChebyshevOperator  spec(this->mspSetup->specSize());

         // Rescale results
         this->mTmpIn.topRows(this->mspSetup->specSize()) = spec.diff(0,1)*chebVal.topRows(this->mspSetup->specSize());

      // Compute simple projection
      } else
      {
         // Copy into other array
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize());
      }

      // Set the padded values to zero
      this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), rPhysVal.data());
   }

   template <Arithmetics::Id TOperation> void ChebyshevFftwTransform::integrate(MatrixZ& rChebVal, const MatrixZ& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rChebVal.rows() == this->mspSetup->bwdSize());
      assert(rChebVal.cols() == this->mspSetup->howmany());

      // Do transform of real part
      this->mTmpIn = physVal.real();
      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());

      // Rescale FFT output
      rChebVal.real() = this->mspSetup->scale()*this->mTmpOut;

      // Do transform of imaginary part
      this->mTmpIn = physVal.imag();
      fftw_execute_r2r(this->mFPlan, this->mTmpIn.data(), this->mTmpOut.data());

      // Rescale FFT output
      rChebVal.imag() = this->mspSetup->scale()*this->mTmpOut;
   }

   template <Arithmetics::Id TOperation> void ChebyshevFftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ChebyshevFftwTransform::ProjectorType::Id projector)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

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
      if(projector == ChebyshevFftwTransform::ProjectorType::DIFF)
      {
         Spectral::ChebyshevOperator  spec(this->mspSetup->specSize());

         // Compute derivative
         this->mTmpIn.topRows(this->mspSetup->specSize()) = spec.diff(0,1)*chebVal.topRows(this->mspSetup->specSize()).real();

      // Compute simple projection of real part
      } else
      {
         // Copy values into simple matrix
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).real();
      }

      // Set the padded values to zero
      this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform of real part
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
      rPhysVal.real() = this->mTmpOut;

      // Compute first derivative of imaginary part
      if(projector == ChebyshevFftwTransform::ProjectorType::DIFF)
      {
         Spectral::ChebyshevOperator  spec(this->mspSetup->specSize());

         // Rescale results
         this->mTmpIn.topRows(this->mspSetup->specSize()) = spec.diff(0,1)*chebVal.topRows(this->mspSetup->specSize()).imag();

      // Compute simple projection of imaginary part
      } else
      {
         // Rescale results
         this->mTmpIn.topRows(this->mspSetup->specSize()) = chebVal.topRows(this->mspSetup->specSize()).imag();
      }

      // Set the padded values to zero
      this->mTmpIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform of imaginary part
      fftw_execute_r2r(this->mBPlan, this->mTmpIn.data(), this->mTmpOut.data());
      rPhysVal.imag() = this->mTmpOut;
   }

}
}

#endif // CHEBYSHEVFFTWTRANSFORM_HPP
