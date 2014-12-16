/**
 * @file ChebyshevFftwTransform.hpp
 * @brief Implementation of the FFTW transform for a Chebyshev expansion 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include <set>
#include <map>

// External includes
//
#include <fftw3.h>

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Exceptions/Exception.hpp"
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
         /** Enum of projector IDs:
          *    - PROJ: projection
          *    - DIFF: D
          */
         enum Id {PROJ, DIFF};
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
          * @brief FFTW plan for the forward transform (real -> real)
          */
         fftw_plan   mFPlan;

         /**
          * @brief FFTW plan for the backward transform (real -> real)
          */
         fftw_plan   mBPlan;

         /**
          * @brief Scale factor for Chebyshev axis
          */
         MHDFloat mCScale;

         /**
          * @brief Storage for data input
          */
         Matrix   mTmpIn;

         /**
          * @brief Storage for data output
          */
         Matrix   mTmpOut;

         #if defined GEOMHDISCC_TRANSOP_FORWARD || defined GEOMHDISCC_TRANSOP_BACKWARD
         /**
          * @brief Storage for the Chebyshev differentiation matrix
          */
         SparseMatrix   mDiff;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD || defined GEOMHDISCC_TRANSOP_BACKWARD

         #if defined GEOMHDISCC_TRANSOP_BACKWARD
         /**
          * @brief Storage for the sparse solver for differentiation
          */
         Solver::SparseSelector<SparseMatrix>::Type mSDiff;

         /**
          * @brief Storage for the backward operators input data
          */
         Matrix mTmpInS;

         /**
          * @brief Storage for the backward operators output data
          */
         Matrix mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_BACKWARD

         /**
          * @brief Initialise the FFTW transforms (i.e. create plans, etc)
          */
         void initFft();

         #if defined GEOMHDISCC_TRANSOP_FORWARD || defined GEOMHDISCC_TRANSOP_BACKWARD
         /**
          * @brief Initialise the spectral operators
          */
         void initOperators();

         #elif defined GEOMHDISCC_TRANSOP_RECURRENCE

         /**
          * @brief Compute derivative by recurrence relation
          */
         void recurrenceDiff(Matrix& rDealiased, const Matrix& rChebVal) const;

         /**
          * @brief Compute derivative by recurrence relation
          */
         void recurrenceDiff(Matrix& rDealiased, const MatrixZ& rChebVal, const bool useImag) const;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD || defined GEOMHDISCC_TRANSOP_BACKWARD


         /**
          * @brief Cleanup memory used by FFTW on destruction
          */
         void cleanupFft();
   };

   template <Arithmetics::Id TOperation> void ChebyshevFftwTransform::integrate(Matrix& rChebVal, const Matrix& physVal, ChebyshevFftwTransform::IntegratorType::Id integrator)
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

      // Do transform
      fftw_execute_r2r(this->mFPlan, const_cast<MHDFloat *>(physVal.data()), rChebVal.data());

      // Rescale to remove FFT scaling
      rChebVal *= this->mspSetup->scale();
   }

   template <Arithmetics::Id TOperation> void ChebyshevFftwTransform::project(Matrix& rPhysVal, const Matrix& chebVal, ChebyshevFftwTransform::ProjectorType::Id projector)
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
      if(projector == ChebyshevFftwTransform::ProjectorType::DIFF)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff*chebVal.topRows(this->mspSetup->specSize());
         #elif defined GEOMHDISCC_TRANSOP_RECURRENCE
            this->recurrenceDiff(this->mTmpIn, chebVal);
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiff, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

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

      // Assert that a mixed transform was setup
      assert(this->mspSetup->type() == FftSetup::COMPONENT);

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
      if(projector == ChebyshevFftwTransform::ProjectorType::DIFF)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).real();
         #elif defined GEOMHDISCC_TRANSOP_RECURRENCE
            this->recurrenceDiff(this->mTmpIn, chebVal, false);
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).real(); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiff, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

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
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiff*chebVal.topRows(this->mspSetup->specSize()).imag();
         #elif defined GEOMHDISCC_TRANSOP_RECURRENCE
            this->recurrenceDiff(this->mTmpIn, chebVal, true);
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).imag(); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiff, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Compute simple projection of imaginary part
      } else
      {
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
