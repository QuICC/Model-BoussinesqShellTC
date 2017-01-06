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

#ifdef QUICC_DEBUG
   #include <helper_cuda.h>
#else
   #define checkCudaErrors(x) x
#endif //QUICC_DEBUG

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/Dimensions.hpp"
#include "Enums/NonDimensional.hpp"
#include "FastTransforms/FftSetup.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "SparseSolvers/SparseLinearSolverTools.hpp"

namespace QuICC {

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
          */
         void integrate(Matrix& rChebVal, const Matrix& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute forward FFT (C2C)
          *
          * Compute the FFT from real physical space to Chebyshev spectral space
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          */
         void integrate(MatrixZ& rChebVal, const MatrixZ& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute backward FFT (R2R)
          *
          * Compute the FFT from Chebyshev spectral space to real physical space
          *
          * @param rPhysVal   Output physical values
          * @param chebVal    Input Chebyshev coefficients
          * @param projector  Projector to use
          */
         void project(Matrix& rPhysVal, const Matrix& chebVal, ProjectorType::Id projector);

         /**
          * @brief Compute backward FFT (C2C)
          *
          * Compute the FFT from Chebyshev spectral space to real physical space
          *
          * @param rPhysVal   Output physical values
          * @param chebVal    Input Chebyshev coefficients
          * @param projector  Projector to use
          */
         void project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ProjectorType::Id projector);

     #ifdef QUICC_STORAGEPROFILE
         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
     #endif // QUICC_STORAGEPROFILE
         
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

}
}

#endif // CHEBYSHEVCUFFTTRANSFORM_HPP
