/**
 * @file FftwTransform.hpp
 * @brief Implementation of the FFTW transform 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FFTWTRANSFORM_HPP
#define FFTWTRANSFORM_HPP

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
#include <fftw3.h>

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
   class FftwTransform 
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
         FftwTransform();

         /**
          * @brief Destroy the FFTW plans
          */
         ~FftwTransform();

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
         fftw_plan   mFPlan;

         /**
          * @brief Plan for the backward transform (complex->real or complex->complex)
          */
         fftw_plan   mBPlan;

         /**
          * @brief Temporary storage used in the projections (complex -> real)
          */
         MatrixZ  mTmpRIn;

         /**
          * @brief Temporary storage used in the projections (complex -> complex)
          */
         MatrixZ  mTmpZIn;

         /**
          * @brief Initialise the FFTW transforms (i.e. create plans, etc)
          */
         void initFft();

         /**
          * @brief Cleanup memory used by FFTW on destruction
          */
         void cleanupFft();
   };

   template <Arithmetics::Id TOperation> void FftwTransform::integrate(MatrixZ& rFFTVal, const Matrix& physVal, FftwTransform::IntegratorType::Id integrator)
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
      fftw_execute_dft_r2c(this->mFPlan, const_cast<MHDFloat *>(physVal.data()), reinterpret_cast<fftw_complex* >(rFFTVal.data()));

      // Rescale output from FFT
      rFFTVal *= this->mspSetup->scale();
   }

   template <Arithmetics::Id TOperation> void FftwTransform::project(Matrix& rPhysVal, const MatrixZ& fftVal, FftwTransform::ProjectorType::Id projector)
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
      if(projector == FftwTransform::ProjectorType::DIFF)
      {
         // Get differentiation factors
         ArrayZ factor = MathConstants::cI*Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);

         // Rescale results
         this->mTmpRIn.topRows(this->mspSetup->specSize()) = factor.asDiagonal()*fftVal.topRows(this->mspSetup->specSize());

      // Compute simple projection
      } else
      {
         // Rescale results
         this->mTmpRIn.topRows(this->mspSetup->specSize()) = fftVal.topRows(this->mspSetup->specSize());
      }

      // Set the padded values to zero
      this->mTmpRIn.row(0).imag().setConstant(0);

      // Set the padded values to zero
      this->mTmpRIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform
      fftw_execute_dft_c2r(this->mBPlan, reinterpret_cast<fftw_complex* >(this->mTmpRIn.data()), rPhysVal.data());
   }

   template <Arithmetics::Id TOperation> void FftwTransform::integrate(MatrixZ& rFFTVal, const MatrixZ& physVal, FftwTransform::IntegratorType::Id integrator)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // Assert that a non mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::EQUAL);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rFFTVal.rows() == this->mspSetup->bwdSize());
      assert(rFFTVal.cols() == this->mspSetup->howmany());

      // Do transform
      fftw_execute_dft(this->mFPlan, reinterpret_cast<fftw_complex* >(const_cast<MHDComplex *>(physVal.data())), reinterpret_cast<fftw_complex* >(rFFTVal.data()));

      // Rescale output from FFT
      rFFTVal *= this->mspSetup->scale();
   }

   template <Arithmetics::Id TOperation> void FftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& fftVal, FftwTransform::ProjectorType::Id projector)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // Assert that a non mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::EQUAL);

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
      if(projector == FftwTransform::ProjectorType::DIFF)
      {
         // Get differentiation factors
         ArrayZ factor = MathConstants::cI*Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);

         // Rescale results
         this->mTmpZIn.topRows(this->mspSetup->specSize()) = factor.asDiagonal()*fftVal.topRows(this->mspSetup->specSize());

      // Compute simple projection
      } else
      {
         // Rescale results
         this->mTmpZIn.topRows(this->mspSetup->specSize()) = fftVal.topRows(this->mspSetup->specSize());
      }

      // Set the padded values to zero
      this->mTmpZIn.bottomRows(this->mspSetup->padSize()).setZero();

      // Do transform
      fftw_execute_dft(this->mBPlan, reinterpret_cast<fftw_complex *>(this->mTmpZIn.data()), reinterpret_cast<fftw_complex *>(rPhysVal.data()));
   }

}
}

#endif // FFTWTRANSFORM_HPP
