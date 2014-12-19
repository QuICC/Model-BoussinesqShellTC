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
          *    - DIFF: D
          *    - DIVR: 1/r
          *    - DIVR2: 1/r^2
          *    - DIVRDIFFR: 1/r D r
          *    - DIFFDIVR: D 1/r
          */
         enum Id {PROJ, DIFF, DIVR, DIVR2, DIVRDIFFR, DIFFDIVR};
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
          * @brief FFTW plan for the even forward transform (real -> real)
          */
         fftw_plan   mFEPlan;

         /**
          * @brief FFTW plan for the odd forward and backward transform (real -> real) 
          *
          * DCT IV is it's own inverse
          */
         fftw_plan   mFBOPlan;

         /**
          * @brief FFTW plan for the even backward transform (real -> real)
          */
         fftw_plan   mBEPlan;

         /**
          * @brief FFTW plan for the even backward transform with odd size (real -> real)
          */
         fftw_plan   mBEOPlan;

         /**
          * @brief FFTW plan for the odd backward transform with even size (real -> real)
          */
         fftw_plan   mFBOEPlan;

         /**
          * @brief Storage for even data input
          */
         Matrix   mTmpEIn;

         /**
          * @brief Storage for odd data input
          */
         Matrix   mTmpOIn;

         /**
          * @brief Storage for even data output
          */
         Matrix   mTmpEOut;

         /**
          * @brief Storage for odd data output
          */
         Matrix   mTmpOOut;

         /**
          * @brief Storage for the even Chebyshev differentiation matrix
          */
         SparseMatrix   mDiffE;

         /**
          * @brief Storage for the odd Chebyshev differentiation matrix
          */
         SparseMatrix   mDiffO;

         #if defined GEOMHDISCC_TRANSOP_FORWARD
         /**
          * @brief Storage for the division by R physical array
          */
         Array   mDivR;

         /**
          * @brief Storage for the division by R^2 physical array
          */
         Array   mDivR2;

         #elif defined GEOMHDISCC_TRANSOP_BACKWARD
         /**
          * @brief Storage for the even division by R matrix
          */
         SparseMatrix   mDivRE;

         /**
          * @brief Storage for the odd division by R matrix
          */
         SparseMatrix   mDivRO;

         /**
          * @brief Storage for the even division by R^2 matrix
          */
         SparseMatrix   mDivR2E;

         /**
          * @brief Storage for the odd division by R^2 matrix
          */
         SparseMatrix   mDivR2O;

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
          * @brief Storage for the sparse solver for even division by R^2
          */
         Solver::SparseSelector<SparseMatrix>::Type mSDivR2E;

         /**
          * @brief Storage for the sparse solver for odd division by R^2
          */
         Solver::SparseSelector<SparseMatrix>::Type mSDivR2O;

         /**
          * @brief Storage for the backward operators input data
          */
         Matrix mTmpInS;

         /**
          * @brief Storage for the backward operators output data
          */
         Matrix mTmpOutS;
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD

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

         /**
          * @brief Extract the parity modes from the whole data
          */
         void extractParityModes(Matrix& rSelected, const Matrix& data, const MatrixI& info, const int rows);

         /**
          * @brief Extract the parity modes from the whole data
          */
         void extractParityModes(Matrix& rSelected, const MatrixZ& data, const bool isReal, const MatrixI& info, const int rows);

         /**
          * @brief Set the parity modes into the whole data
          */
         void setParityModes(Matrix& rData, const Matrix& selected, const MatrixI& info, const int rows);

         /**
          * @brief Set the parity modes into the whole data
          */
         void setParityModes(MatrixZ& rData, const Matrix& selected, const bool isReal, const MatrixI& info, const int rows);

         /**
          * @brief Add the parity modes into the whole data
          */
         void addParityModes(Matrix& rData, const Matrix& selected, const MatrixI& info, const int rows);

         /**
          * @brief Add the parity modes into the whole data
          */
         void addParityModes(MatrixZ& rData, const Matrix& selected, const bool isReal, const MatrixI& info, const int rows);

         /**
          * @brief Get input temporary data depending on parity
          */
         Matrix& rTmpIn(const int parity);

         /**
          * @brief Get output temporary data depending on parity
          */
         Matrix& rTmpOut(const int parity);

         /**
          * @brief Get the parity block description
          */
         const MatrixI& parityBlocks(const int parity) const;

         /**
          * @brief Get the FFTW plan depending on parity
          */
         fftw_plan fPlan(const int parity);

         /**
          * @brief Get the FFTW plan depending on parity
          */
         fftw_plan bPlan(const int parity, const int sizeParity);

         /**
          * @brief Get the differentiation matrix depending on parity
          */
         const SparseMatrix& diff(const int parity) const;
   };

   inline Matrix& CylinderChebyshevFftwTransform::rTmpIn(const int parity)
   {
      if(parity == 0)
      {
         return this->mTmpEIn;
      } else
      {
         return this->mTmpOIn;
      }
   }

   inline Matrix& CylinderChebyshevFftwTransform::rTmpOut(const int parity)
   {
      if(parity == 0)
      {
         return this->mTmpEOut;
      } else
      {
         return this->mTmpOOut;
      }
   }

   inline const MatrixI& CylinderChebyshevFftwTransform::parityBlocks(const int parity) const
   {
      if(parity == 0)
      {
         return this->mspSetup->evenBlocks();
      } else
      {
         return this->mspSetup->oddBlocks();
      }
   }

   inline fftw_plan CylinderChebyshevFftwTransform::fPlan(const int parity)
   {
      if(parity == 0)
      {
         return this->mFEPlan;
      } else
      {
         return this->mFBOPlan;
      }
   }

   inline fftw_plan CylinderChebyshevFftwTransform::bPlan(const int parity, const int sizeParity)
   {
      if(parity == 0 && parity == sizeParity)
      {
         return this->mBEPlan;
      } else if(parity == 1 && parity == sizeParity)
      {
         return this->mFBOPlan;
      } else if(parity == 0)
      {
         return this->mBEOPlan;
      } else
      {
         return this->mFBOEPlan;
      }
   }

   inline const SparseMatrix& CylinderChebyshevFftwTransform::diff(const int parity) const
   {
      if(parity == 0)
      {
         return this->mDiffE;
      } else
      {
         return this->mDiffO;
      }
   }

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

      // Do even transform
      for(int parity = 0; parity < 2; ++parity)
      {
         this->extractParityModes(this->rTmpIn(parity), physVal, this->parityBlocks(parity), physVal.rows());
         fftw_execute_r2r(this->fPlan(parity), this->rTmpIn(parity).data(), this->rTmpOut(parity).data());
         this->setParityModes(rChebVal, this->rTmpOut(parity), this->parityBlocks(parity), physVal.rows());
      }

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

      for(int parity = 0; parity < 2; ++parity)
      {
         int fftParity = parity;

         this->extractParityModes(this->rTmpIn(parity), chebVal, this->parityBlocks(parity), this->mspSetup->specSize());

         // Compute first derivative
         if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIFF)
         {
            fftParity = (fftParity + 1) % 2;

            #if defined GEOMHDISCC_TRANSOP_FORWARD
               this->rTmpIn(parity).topRows(this->mspSetup->specSize()) = this->diff(parity)*this->rTmpIn(parity).topRows(this->mspSetup->specSize());
            #elif defined GEOMHDISCC_TRANSOP_BACKWARD
               throw Exception("Not yet implemented!");
            #endif //defined GEOMHDISCC_TRANSOP_FORWARD

         #if defined GEOMHDISCC_TRANSOP_BACKWARD
         // Compute division by R
         } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR)
         {
            fftParity = (fftParity + 1) % 2;
            throw Exception("Not yet implemented!");

         // Compute division by R^2
         } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR2)
         {
            throw Exception("Not yet implemented!");
         #endif //defined GEOMHDISCC_TRANSOP_BACKWARD

         // Compute 1/r D r
         } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
         {
            throw Exception("DIVRDIFFR operator is not yet implemented");

         // Compute D 1/r
         } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIFFDIVR)
         {
            fftParity = (fftParity + 1) % 2;

            #if defined GEOMHDISCC_TRANSOP_FORWARD
               // Compute D f/r part
               this->rTmpIn(parity).topRows(this->mspSetup->specSize()) = this->diff(parity)*this->rTmpIn(parity).topRows(this->mspSetup->specSize());
            #elif defined GEOMHDISCC_TRANSOP_BACKWARD
               throw Exception("Backward DIFFDIVR operator is not yet implemented");
            #endif //defined GEOMHDISCC_TRANSOP_FORWARD
         }

         // Set the padded values to zero
         this->rTmpIn(parity).bottomRows(this->mspSetup->padSize()).setZero();

         // Do transform
         fftw_execute_r2r(this->bPlan(fftParity,parity), this->rTmpIn(parity).data(), this->rTmpOut(parity).data());

         #if defined GEOMHDISCC_TRANSOP_FORWARD
         // Compute division by R
         if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR)
         {
            this->rTmpOut(parity) = this->mDivR.asDiagonal()*this->rTmpOut(parity);

            this->setParityModes(rPhysVal, this->rTmpOut(parity), this->parityBlocks(parity), rPhysVal.rows());

         // Compute division by R^2
         } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR2)
         {
            this->rTmpOut(parity) = this->mDivR2.asDiagonal()*this->rTmpOut(parity);

            this->setParityModes(rPhysVal, this->rTmpOut(parity), this->parityBlocks(parity), rPhysVal.rows());

         // Compute D 1/r
         } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIFFDIVR)
         {
            // Compute Df/r
            this->rTmpOut(parity) = this->mDivR.asDiagonal()*this->rTmpOut(parity);

            this->setParityModes(rPhysVal, this->rTmpOut(parity), this->parityBlocks(parity), rPhysVal.rows());

            // Compute f/r^2
            this->extractParityModes(this->rTmpIn(parity), chebVal, this->parityBlocks(parity), this->mspSetup->specSize());

            // Set the padded values to zero
            this->rTmpIn(parity).bottomRows(this->mspSetup->padSize()).setZero();

            // Do transform
            fftw_execute_r2r(this->bPlan(parity,parity), this->rTmpIn(parity).data(), this->rTmpOut(parity).data());

            // Compute D f/r - f/r^2
            this->rTmpOut(parity) = (-this->mDivR2).asDiagonal()*this->rTmpOut(parity);

            this->addParityModes(rPhysVal, this->rTmpOut(parity), this->parityBlocks(parity), rPhysVal.rows());
         } else
         {
            this->setParityModes(rPhysVal, this->rTmpOut(parity), this->parityBlocks(parity), rPhysVal.rows());
         }
         #else
            this->setParityModes(rPhysVal, this->rTmpOut(parity), this->parityBlocks(parity), rPhysVal.rows());
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD
      }
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

      // Loop over real and imaginary
      for(int isReal = 0; isReal < 2; ++isReal)
      {
         // Do transform of real part
         for(int parity = 0; parity < 2; ++parity)
         {
            this->extractParityModes(this->rTmpIn(parity), physVal, isReal, this->parityBlocks(parity), physVal.rows());
            fftw_execute_r2r(this->fPlan(parity), this->rTmpIn(parity).data(), this->rTmpOut(parity).data());
            this->setParityModes(rChebVal, this->rTmpOut(parity), isReal, this->parityBlocks(parity), physVal.rows());
         }
      }

      // Rescale FFT output
      rChebVal *= this->mspSetup->scale();
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

      // Loop over real and imaginary
      for(int isReal = 0; isReal < 2; ++isReal)
      {
         for(int parity = 0; parity < 2; ++parity)
         {
            int fftParity = parity;

            this->extractParityModes(this->rTmpIn(parity), chebVal, isReal, this->parityBlocks(parity), this->mspSetup->specSize());

            // Compute first derivative of real part
            if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIFF)
            {
               fftParity = (fftParity + 1) % 2;

               #if defined GEOMHDISCC_TRANSOP_FORWARD
                  this->rTmpIn(parity).topRows(this->mspSetup->specSize()) = this->diff(parity)*this->rTmpIn(parity).topRows(this->mspSetup->specSize());
               #elif defined GEOMHDISCC_TRANSOP_BACKWARD
                  throw Exception("Not yet implemented!");
               #endif //defined GEOMHDISCC_TRANSOP_FORWARD

            #if defined GEOMHDISCC_TRANSOP_BACKWARD
            // Compute division by R of real part
            } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR)
            {
               fftParity = (fftParity + 1) % 2;
               throw Exception("Not yet implemented!");

            // Compute division by R^2 of real part
            } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR2)
            {
               throw Exception("Not yet implemented!");
            #endif //defined GEOMHDISCC_TRANSOP_BACKWARD

            // Compute 1/r D r
            } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVRDIFFR)
            {
               throw Exception("DIVRDIFFR operator is not yet implemented");

            // Compute D 1/r
            } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIFFDIVR)
            {
               fftParity = (fftParity + 1) % 2;

               #if defined GEOMHDISCC_TRANSOP_FORWARD
                  // Compute D f/r part
                  this->rTmpIn(parity).topRows(this->mspSetup->specSize()) = this->diff(parity)*this->rTmpIn(parity).topRows(this->mspSetup->specSize());
               #elif defined GEOMHDISCC_TRANSOP_BACKWARD
                  throw Exception("Backward DIFFDIVR operator is not yet implemented");
               #endif //defined GEOMHDISCC_TRANSOP_FORWARD
            }

            // Set the padded values to zero
            this->rTmpIn(parity).bottomRows(this->mspSetup->padSize()).setZero();

            // Do transform of real part
            fftw_execute_r2r(this->bPlan(fftParity,parity), this->rTmpIn(parity).data(), this->rTmpOut(parity).data());

            #if defined GEOMHDISCC_TRANSOP_FORWARD
            // Compute division by R
            if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR)
            {
               this->rTmpOut(parity) = this->mDivR.asDiagonal()*this->rTmpOut(parity);

               this->setParityModes(rPhysVal, this->rTmpOut(parity), isReal, this->parityBlocks(parity), rPhysVal.rows());

            // Compute division by R^2
            } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIVR2)
            {
               this->rTmpOut(parity) = this->mDivR2.asDiagonal()*this->rTmpOut(parity);

               this->setParityModes(rPhysVal, this->rTmpOut(parity), isReal, this->parityBlocks(parity), rPhysVal.rows());

            // Compute D 1/r
            } else if(projector == CylinderChebyshevFftwTransform::ProjectorType::DIFFDIVR)
            {
               // Compute Df/r
               this->rTmpOut(parity) = this->mDivR.asDiagonal()*this->rTmpOut(parity);

               this->setParityModes(rPhysVal, this->rTmpOut(parity), isReal, this->parityBlocks(parity), rPhysVal.rows());

               // Compute f/r^2
               this->extractParityModes(this->rTmpIn(parity), chebVal, isReal, this->parityBlocks(parity), this->mspSetup->specSize());

               // Set the padded values to zero
               this->rTmpIn(parity).bottomRows(this->mspSetup->padSize()).setZero();

               // Do transform
               fftw_execute_r2r(this->bPlan(parity,parity), this->rTmpIn(parity).data(), this->rTmpOut(parity).data());
               // Compute D f/r - f/r^2
               this->rTmpOut(parity) = (-this->mDivR2).asDiagonal()*this->rTmpOut(parity);

               this->addParityModes(rPhysVal, this->rTmpOut(parity), isReal, this->parityBlocks(parity), rPhysVal.rows());
            } else
            {
               this->setParityModes(rPhysVal, this->rTmpOut(parity), isReal, this->parityBlocks(parity), rPhysVal.rows());
            }
            #else
               this->setParityModes(rPhysVal, this->rTmpOut(parity), isReal, this->parityBlocks(parity), rPhysVal.rows());
         #endif //defined GEOMHDISCC_TRANSOP_FORWARD
         }
      }
   }

}
}

#endif // CYLINDERCHEBYSHEVFFTWTRANSFORM_HPP
