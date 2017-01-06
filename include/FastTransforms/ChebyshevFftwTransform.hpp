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
         /** Enum of projector IDs:
          *    - PROJ: projection
          *    - DIFF: D
          *    - DIFF2: D^2
          */
         enum Id {PROJ, DIFF, DIFF2, DIFFM};
      };

      /**
       * @brief Simple struct holding the integrator IDs
       */
      struct Integrators
      {
         /// Enum of integrator IDs
         // INTG: Standard Chebyshev integrator
         // INTGI2: Chebyshev integrator with i1
         // INTGI4: Chebyshev integrator with i4
         // INTGI2D1: Chebyshev integrator with i2d1
         // INTGI4D1: Chebyshev integrator with i4d1
         // INTGI2D1MI2: Chebyshev integrator with i2d1 and i2 for mean
         // INTGI4D1MI2: Chebyshev integrator with i4d1 and i2 for mean
         enum Id {INTG, INTGI2, INTGI4, INTGI2D1, INTGI4D1, INTGI2D1MI2, INTGI4D1MI2};
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

         /**
          * @brief Compute forward FFT (R2R) provide full output without spectral truncation
          *
          * Compute the FFT from real physical space to Chebyshev spectral space
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          */
         void integrate_full(Matrix& rChebVal, const Matrix& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute forward FFT (C2C) provide full output without spectral truncation
          *
          * Compute the FFT from real physical space to Chebyshev spectral space
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          */
         void integrate_full(MatrixZ& rChebVal, const MatrixZ& physVal, IntegratorType::Id integrator);

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

         /**
          * @brief Storage for the projector operators
          */
         std::map<ProjectorType::Id, SparseMatrix> mProjOp;

         /**
          * @brief Storage for the integrator operators
          */
         std::map<IntegratorType::Id, SparseMatrix> mIntgOp;

         /**
          * @brief Storage for the sparse solver matrices
          */
         std::map<ProjectorType::Id, SparseMatrix> mSolveOp;

         /**
          * @brief Storage for the sparse triangular solvers
          */
         std::map<ProjectorType::Id, SharedPtrMacro<Solver::SparseTriSelector<SparseMatrix>::Type> > mTriSolver;

         /**
          * @brief Storage for the sparse SPD solvers
          */
         std::map<ProjectorType::Id, SharedPtrMacro<Solver::SparseSpdSelector<SparseMatrix>::Type> > mSpdSolver;

         /**
          * @brief Storage for the backward operators input data
          */
         Matrix mTmpInS;

         /**
          * @brief Storage for the backward operators output data
          */
         Matrix mTmpOutS;

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

}
}

#endif // CHEBYSHEVFFTWTRANSFORM_HPP
