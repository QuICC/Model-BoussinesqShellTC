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
#include "Enums/Dimensions.hpp"
#include "Enums/NonDimensional.hpp"
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
         // PROJ: Projector
         // DIFF: First derivative
         // DIFF2: Second derivative
         // DIFF3: Third derivative
         // PROJMEANONLY: Only project mean component (WARNING: implementation is slow)
         enum Id {PROJ, DIFF, DIFF2, DIFF3, PROJMEANONLY};
      };

      /**
       * @brief Simple struct holding the integrator IDs
       */
      struct Integrators
      {
         /// Enum of integrator IDs
         // INTG: Integrator
         // INTGM: Integrator and zero k2 = 0, k1 != 0 modes (only complex - complex)
         // INTGDIFF: First derivative integrator
         // INTGDIFF2: Second derivative integrator
         // INTGDIFFM: First derivative integrator and mean (k = 0 mode is not zeroed)
         // INTGDIFFNEGM: First derivative integrator and negative mean (k = 0 mode is not zeroed)
         // INTGINVLAPLH: Inverse horizontal laplacian
         // INTGINVDIFFLAPLHH: Derivative of inverse horizontal laplacian
         // INTGMEANONLY: Only compute the mean component (WARNING implementation is slow!)
         enum Id {INTG, INTGDIFF, INTGDIFF2, INTGM, INTGDIFFM, INTGDIFFNEGM, INTGINVLAPLH, INTGDIFFINVLAPLH, INTGMEANONLY};
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

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief Impose complex conjugate symmetry on special indexes
          */
         void forceConjugate(MatrixZ& rFFTVal);

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
          * @brief Temporary storage used in the projections (complex -> real)
          */
         Matrix  mTmpROut;

         /**
          * @brief Temporary storage used in the projections (complex -> complex)
          */
         MatrixZ  mTmpZIn;

         /**
          * @brief Temporary storage used in the projections (complex -> complex)
          */
         MatrixZ  mTmpZOut;

         /**
          * @brief Storage for the mean block sizes
          */
         std::vector<std::pair<int,int> > mMeanBlocks;

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

#endif // FFTWTRANSFORM_HPP
