/**
 * @file CylinderChebyshevFftwTransform.hpp
 * @brief Implementation of the FFTW transform for a Chebyshev expansion for a cylinder radius  
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CYLINDERCHEBYSHEVFFTWTRANSFORM_HPP
#define CYLINDERCHEBYSHEVFFTWTRANSFORM_HPP

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
    * @brief Simple struct holding details about CylinderChebyshevFFT transform
    */
   struct CylinderChebyshevFftIds {

      /**
       * @brief Simple struct holding the projector IDs
       */
      struct Projectors
      {
         /** Enum of projector IDs:
          *    - PROJ: projection
          *    - DIFF: derivative
          *    - DIVR: division by R
          */
         enum Id {PROJ,  DIFF, DIVR};
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
    * @brief Implementation of the FFTW transform for a Chebyshev expansion for a cylinder radius
    */ 
   class CylinderChebyshevFftwTransform
   {
      public:
         /// Typedef for the configuration class
         typedef FftSetup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedFftSetup SharedSetupType;

         /// Typedef for the Projector type
         typedef CylinderChebyshevFftIds::Projectors ProjectorType;

         /// Typedef for the Integrator type
         typedef CylinderChebyshevFftIds::Integrators IntegratorType;

         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size); 

         /**
          * @brief Very basic constructor
          */
         CylinderChebyshevFftwTransform();

         /**
          * @brief Destroy the FFTW plans
          */
         ~CylinderChebyshevFftwTransform();

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
          * \mhdTodo Projectors should be converted to Python
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
          * \mhdTodo Projectors should be converted to Python
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
          * @brief Storage for the even Chebyshev differentiation matrix
          */
         SparseMatrix   mDiffE;

         /**
          * @brief Storage for the odd Chebyshev differentiation matrix
          */
         SparseMatrix   mDiffO;

         /**
          * @brief Storage for the even division by R matrix
          */
         SparseMatrix   mDivRE;

         /**
          * @brief Storage for the odd division by R matrix
          */
         SparseMatrix   mDivRO;

         #if defined GEOMHDISCC_TRANSOP_BACKWARD
         /**
          * @brief Storage for the sparse solver for even differentiation
          */
         Solver::SparseSelector<SparseMatrix>::Type mSDiffE;

         /**
          * @brief Storage for the sparse solver for odd differentiation
          */
         Solver::SparseSelector<SparseMatrix>::Type mSDiffO;

         /**
          * @brief Storage for the sparse solver for even division by R
          */
         Solver::SparseSelector<SparseMatrix>::Type mSDivRE;

         /**
          * @brief Storage for the sparse solver for odd division by R
          */
         Solver::SparseSelector<SparseMatrix>::Type mSDivRO;

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

         /**
          * @brief Initialise the spectral operators
          */
         void initOperators();

         /**
          * @brief Cleanup memory used by FFTW on destruction
          */
         void cleanupFft();
   };

   template <Arithmetics::Id TOperation> void CylinderChebyshevFftwTransform::integrate(Matrix& rChebVal, const Matrix& physVal, CylinderChebyshevFftwTransform::IntegratorType::Id integrator)
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

   template <Arithmetics::Id TOperation> void CylinderChebyshevFftwTransform::project(Matrix& rPhysVal, const Matrix& chebVal, CylinderChebyshevFftwTransform::ProjectorType::Id projector)
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
      if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIFF)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiffO*chebVal.topRows(this->mspSetup->specSize());
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiffO, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Compute division by R
      } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDivRO*chebVal.topRows(this->mspSetup->specSize());
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()); 
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDivRO, this->mTmpInS);
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

   template <Arithmetics::Id TOperation> void CylinderChebyshevFftwTransform::integrate(MatrixZ& rChebVal, const MatrixZ& physVal, CylinderChebyshevFftwTransform::IntegratorType::Id integrator)
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

   template <Arithmetics::Id TOperation> void CylinderChebyshevFftwTransform::project(MatrixZ& rPhysVal, const MatrixZ& chebVal, CylinderChebyshevFftwTransform::ProjectorType::Id projector)
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
      if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIFF)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiffO*chebVal.topRows(this->mspSetup->specSize()).real();
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).real(); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiffO, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Compute division by R of real part
      } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDivRO*chebVal.topRows(this->mspSetup->specSize()).real();
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).real(); 
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDivRO, this->mTmpInS);
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
      if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIFF)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDiffO*chebVal.topRows(this->mspSetup->specSize()).imag();
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).imag(); 
            this->mTmpInS.topRows(1).setZero();
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDiffO, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

      // Compute division by R of imaginary part
      } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR)
      {
         #if defined GEOMHDISCC_TRANSOP_FORWARD
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mDivRO*chebVal.topRows(this->mspSetup->specSize()).imag();
         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
            this->mTmpInS = chebVal.topRows(this->mspSetup->specSize()).imag(); 
            Solver::internal::solveWrapper(this->mTmpOutS, this->mSDivRO, this->mTmpInS);
            this->mTmpIn.topRows(this->mspSetup->specSize()) = this->mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

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

#endif // CYLINDERCHEBYSHEVFFTWTRANSFORM_HPP
