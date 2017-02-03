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
#include "FastTransforms/CuFftLibrary.hpp"

namespace QuICC {

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
          * @brief Compute forward FFT (R2C)
          *
          * Compute the FFT from real physical space to complex spectral space
          *
          * @param rFFTVal Output FFT transformed values
          * @param physVal Input physical values
          * @param integrator Integrator to use
          */
         void integrate(MatrixZ& rFFTVal, const Matrix& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute backward FFT (C2R)
          *
          * Compute the FFT from complex spectral space to real physical space
          *
          * @param rPhysVal Output physical values
          * @param fftVal Input FFT values
          * @param projector  Projector to use
          */
         void project(Matrix& rPhysVal, const MatrixZ& fftVal, ProjectorType::Id projector);

         /**
          * @brief Compute forward FFT (C2C)
          *
          * Compute the FFT from complex physical space to complex spectral space
          *
          * @param rFFTVal Output FFT transformed values
          * @param physVal Input physical values
          * @param integrator Integrator to use
          */
         void integrate(MatrixZ& rFFTVal, const MatrixZ& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute backward FFT (C2C)
          *
          * Compute the FFT from complex spectral space to complex physical space
          *
          * @param rPhysVal   Output physical values
          * @param fftVal     Input FFT values
          * @param projector  Projector to use
          */
         void project(MatrixZ& rPhysVal, const MatrixZ& fftVal, ProjectorType::Id projector);

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
          * @brief Batch size for each CUDA stream
          */
         std::vector<int>  mStreamBatch;

         /**
          * @brief Plan for the forward transform (real->complex or complex->complex)
          */
         std::vector<cufftHandle>   mFPlan;

         /**
          * @brief Plan for the backward transform (complex->real or complex->complex)
          */
         std::vector<cufftHandle>   mBPlan;

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
         std::vector<cufftDoubleReal *> mpDevR;

         /**
          * @brief Complex storage on device
          */
         std::vector<cufftDoubleComplex *> mpDevZI;

         /**
          * @brief Complex storage on device
          */
         std::vector<cufftDoubleComplex *> mpDevZO;

         /**
          * @brief Initialise the FFTW transforms (i.e. create plans, etc)
          */
         void initFft();

         /**
          * @brief Cleanup memory used by FFTW on destruction
          */
         void cleanupFft();
   };

}
}

#endif // CUFFTTRANSFORM_HPP
