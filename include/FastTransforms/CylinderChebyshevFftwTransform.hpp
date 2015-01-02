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
         /** 
          * Enum of projector IDs:
          *    - PROJ: projection
          *    - DIFF: D
          *    - DIVR: 1/r
          *    - DIVR2: 1/r^2
          *    - DIVRDIFFR: 1/r D r
          *    - DIFFDIVR: D 1/r
          *
          * All projectors assume "even" data, ie even "m" -> even "r", odd "m" -> odd "r"
          */
         enum Id {PROJ, DIFF, DIVR, DIVR2, DIVRDIFFR, DIFFDIVR};
      };

      /**
       * @brief Simple struct holding the integrator IDs
       */
      struct Integrators
      {
         /** 
          * Enum of integrator IDs
          *    - INTGE: even integration, even "m" -> even "r", etc
          *    - INTGO: odd integration, even "m" -> odd "r", etc
          */
         enum Id {INTGE, INTGO};
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
          * @param arithId    Arithmetic operation to perform
          */
         void integrate(Matrix& rChebVal, const Matrix& physVal, IntegratorType::Id integrator, Arithmetics::Id arithId);

         /**
          * @brief Compute forward FFT (C2C)
          *
          * Compute the FFT from real physical space to Chebyshev spectral space
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          * @param arithId    Arithmetic operation to perform
          */
         void integrate(MatrixZ& rChebVal, const MatrixZ& physVal, IntegratorType::Id integrator, Arithmetics::Id arithId);

         /**
          * @brief Compute backward FFT (R2R)
          *
          * Compute the FFT from Chebyshev spectral space to real physical space
          *
          * @param rPhysVal   Output physical values
          * @param chebVal    Input Chebyshev coefficients
          * @param projector  Projector to use
          * @param arithId    Arithmetic operation to perform
          */
         void project(Matrix& rPhysVal, const Matrix& chebVal, ProjectorType::Id projector, Arithmetics::Id arithId);

         /**
          * @brief Compute backward FFT (C2C)
          *
          * Compute the FFT from Chebyshev spectral space to real physical space
          *
          * @param rPhysVal   Output physical values
          * @param chebVal    Input Chebyshev coefficients
          * @param projector  Projector to use
          * @param arithId    Arithmetic operation to perform
          */
         void project(MatrixZ& rPhysVal, const MatrixZ& chebVal, ProjectorType::Id projector, Arithmetics::Id arithId);

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
          * @brief FFTW plan for the odd forward transform with even size (real -> real)
          */
         fftw_plan   mFEOPlan;

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
         fftw_plan fPlan(const int parity, const int sizeParity);

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

   inline fftw_plan CylinderChebyshevFftwTransform::fPlan(const int parity, const int sizeParity)
   {
      if(parity == 0 && parity == sizeParity)
      {
         return this->mFEPlan;
      } else if(parity == 1 && parity == sizeParity)
      {
         return this->mFBOPlan;
      } else if(parity == 0)
      {
         return this->mFEOPlan;
      } else
      {
         return this->mFBOEPlan;
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

}
}

#endif // CYLINDERCHEBYSHEVFFTWTRANSFORM_HPP
