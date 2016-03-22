/**
 * @file SphereChebyshevFftwTransform.hpp
 * @brief Implementation of the FFTW transform for a Chebyshev expansion for a sphere radius 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHERECHEBYSHEVFFTWTRANSFORM_HPP
#define SPHERECHEBYSHEVFFTWTRANSFORM_HPP

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
#include "Enums/Dimensions.hpp"
#include "Enums/Arithmetics.hpp"
#include "Enums/NonDimensional.hpp"
#include "FastTransforms/FftSetup.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "SparseSolvers/SparseLinearSolverTools.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Simple struct holding details about SphercialChebyshevFFT transform
    */
   struct SphereChebyshevFftIds {

      /**
       * @brief Simple struct holding the projector IDs
       */
      struct Projectors
      {
         /** Enum of projector IDs
          *    - PROJ: projection
          *    - DIVR: 1/r
          *    - DIVR2: 1/r^2
          *    - DIFF: D
          *    - DIFF2: D^2
          *    - DIFFR: D r
          *    - DIVRDIFFR: 1/r D r
          *    - RADLAPL: radial laplacian: D^2 + 2/r D
          */
         enum Id {PROJ, DIVR, DIVR2, DIFF, DIFF2, DIFFR, DIVRDIFFR, RADLAPL};
      };

      /**
       * @brief Simple struct holding the integrator IDs
       */
      struct Integrators
      {
         /** 
          * Enum of integrator IDs:
          *    - INTG: integration
          *    - INTGR: integration of r
          *    - INTGQ4: integration of QST Q component for Poloidal NL (4th order equation)
          *    - INTGS4: integration of QST S component for Poloidal NL (4th order equation)
          *    - INTGT: integration of QST T component for Toroidal NL (2nd order equation)
          *    - INTGQ2: integration of QST Q component for Poloidal NL (2nd order equation)
          *    - INTGS2: integration of QST S component for Poloidal NL (2nd order equation)
          */
         enum Id {INTG, INTGR, INTGQ4, INTGS4, INTGT, INTGQ2, INTGS2};
      };

      /**
       * @brief Simple struct holding regularity stencil IDs
       */
      struct Stencils
      {
         /**
          * @brief Enum of regularity constraints:
          *    - REG0:  0th order regularity (i.e. f(0) = 0)
          *    - REG1:  0th order regularity (i.e. f(0) = f'(0) = 0)
          *    - REG2:  0th order regularity (i.e. f(0) = f'(0) = f''(0) = 0)
          *    - REG3:  0th order regularity (i.e. f(0) = f'(0) = f''(0) = f'''(0) = 0)
          *    - REG4:  0th order regularity (i.e. f(0) = f'(0) = f''(0) = f'''(0) = f''''(0) = 0)
          *    - REG5:  0th order regularity (i.e. f(0) = f'(0) = f''(0) = f'''(0) = f''''(0) = ... = 0)
          *    - REG6:  0th order regularity (i.e. f(0) = f'(0) = f''(0) = f'''(0) = f''''(0) = ... = 0)
          *    - REG7:  0th order regularity (i.e. f(0) = f'(0) = f''(0) = f'''(0) = f''''(0) = ... = 0)
          *    - REG8:  0th order regularity (i.e. f(0) = f'(0) = f''(0) = f'''(0) = f''''(0) = ... = 0)
          */
         enum Id {
            REG0 = 0,
            REG1,
            REG2, 
            REG3,
            REG4,
            REG5,
            REG6,
            REG7,
            REG8,
            REGMAX = REG2
         };
      };
   };

   /**
    * @brief Implementation of the FFTW transform for a Chebyshev expansion for a sphere radius
    */ 
   class SphereChebyshevFftwTransform
   {
      public:
         /// Typedef for the configuration class
         typedef FftSetup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedFftSetup SharedSetupType;

         /// Typedef for the Projector type
         typedef SphereChebyshevFftIds::Projectors ProjectorType;

         /// Typedef for the Integrator type
         typedef SphereChebyshevFftIds::Integrators IntegratorType;

         /// Typedef for the Integrator type
         typedef SphereChebyshevFftIds::Stencils RegularityType;

         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size); 

         /**
          * @brief Very basic constructor
          */
         SphereChebyshevFftwTransform();

         /**
          * @brief Destroy the FFTW plans
          */
         ~SphereChebyshevFftwTransform();

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

         /**
          * @brief Compute forward FFT (R2R) provide full output without spectral truncation
          *
          * Compute the FFT from real physical space to Chebyshev spectral space
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          * @param arithId    Arithmetic operation to perform
          */
         void integrate_full(Matrix& rChebVal, const Matrix& physVal, IntegratorType::Id integrator, Arithmetics::Id arithId);

         /**
          * @brief Compute forward FFT (C2C)
          *
          * Compute the FFT from real physical space to Chebyshev spectral space full output without spectral truncation
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          * @param arithId    Arithmetic operation to perform
          */
         void integrate_full(MatrixZ& rChebVal, const MatrixZ& physVal, IntegratorType::Id integrator, Arithmetics::Id arithId);

         /**
          * @brief Compute forward FFT (R2R) for energy calculation on full output without spectral truncation
          *
          * Compute the FFT from real physical space to Chebyshev spectral space for energy (f^2 and thus always even)
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          * @param arithId    Arithmetic operation to perform
          */
         void integrate_energy(Matrix& rChebVal, const Matrix& physVal, IntegratorType::Id integrator, Arithmetics::Id arithId);

         /**
          * @brief Compute forward FFT (C2C) for energy calculation  on full output without spectral truncation
          *
          * Compute the FFT from real physical space to Chebyshev spectral space for energy (f^2 and thus always even)
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          * @param arithId    Arithmetic operation to perform
          */
         void integrate_energy(MatrixZ& rChebVal, const MatrixZ& physVal, IntegratorType::Id integrator, Arithmetics::Id arithId);

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
          * @brief Storage for data input
          */
         Matrix   mTmpIn;

         /**
          * @brief Storage for data output
          */
         Matrix   mTmpOut;

         /**
          * @brief Regularity operators basis size
          */
         std::map<RegularityType::Id, std::pair<int, int> > mRegularitySize;

         /**
          * @brief Does projector flip parity
          */
         std::map<ProjectorType::Id, int> mProjectorFlips;

         /**
          * @brief Does integrator flip parity
          */
         std::map<IntegratorType::Id, int> mIntegratorFlips;

         /**
          * @brief Storage for the regularity constraints operator
          */
         std::map<RegularityType::Id, std::pair<SparseMatrix,SparseMatrix> > mRegOp;

         /**
          * @brief Storage for the sparse solver matrices for the regularity constraints
          */
         std::map<RegularityType::Id, std::pair<SparseMatrix,SparseMatrix> > mRegSolveOp;

         /**
          * @brief Storage for the sparse solvers for the regularity constraints
          */
         std::map<RegularityType::Id, std::pair<SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>,SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type> > > mRegSolver;

         /**
          * @brief Storage for the projector operators
          */
         std::map<ProjectorType::Id, std::pair<SparseMatrix,SparseMatrix> > mProjOp;

         /**
          * @brief Storage for the integrator operators
          */
         std::map<IntegratorType::Id, std::pair<SparseMatrix,SparseMatrix> > mIntgOp;

         /**
          * @brief Storage for the sparse solver matrices
          */
         std::map<ProjectorType::Id, std::pair<SparseMatrix,SparseMatrix> > mSolveOp;

         /**
          * @brief Storage for the sparse solvers
          */
         std::map<ProjectorType::Id, std::pair<SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type>,SharedPtrMacro<Solver::SparseSelector<SparseMatrix>::Type> > > mSolver;

         /**
          * @brief Storage for the backward operators input data for even parity
          */
         Matrix mTmpInSE;

         /**
          * @brief Storage for the backward operators input data for odd parity
          */
         Matrix mTmpInSO;

         /**
          * @brief Storage for the backward operators output data for even parity
          */
         Matrix mTmpOutSE;

         /**
          * @brief Storage for the backward operators output data for odd parity
          */
         Matrix mTmpOutSO;

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
          * @brief Get sparse integration operator
          */
         SparseMatrix& intgOp(const IntegratorType::Id integrator, const int parity);

         /**
          * @brief Get sparse regularity operator
          */
         SparseMatrix& regOp(const RegularityType::Id reg, const int parity);

         /**
          * @brief Get sparse solver operator
          */
         SparseMatrix& solveOp(const ProjectorType::Id projector, const int parity);

         /**
          * @brief Get sparse regularity solver operator
          */
         SparseMatrix& regSolveOp(const RegularityType::Id reg, const int parity);

         /**
          * @brief Get sparse solver
          */
         Solver::SparseSelector<SparseMatrix>::Type& solver(const ProjectorType::Id projector, const int parity);

         /**
          * @brief Get sparse regularity solver
          */
         Solver::SparseSelector<SparseMatrix>::Type& regSolver(const RegularityType::Id reg, const int parity);

         /**
          * @brief Get temporary storage for solver
          */
         Matrix& tmpInS(const int parity);

         /**
          * @brief Get temporary storage for solver
          */
         Matrix& tmpOutS(const int parity);

         /**
          * brief Size of regularity basis
          */
         int regularitySize(const RegularityType::Id reg, const int parity) const;

         /**
          * brief Does operator flip parity?
          */
         int flipsParity(const ProjectorType::Id projector) const;

         /**
          * brief Does operator flip parity?
          */
         int flipsParity(const IntegratorType::Id integrator) const; 

         /**
          * @brief Apply regularity filters on block of data
          */
         void regularizeBlock(Matrix& rData, const int start, const int cols, const int parity, const RegularityType::Id reg);

         /**
          * @brief Apply regularity filters
          */
         void regularize(Matrix& rData, const int parity);
   };

   inline const MatrixI& SphereChebyshevFftwTransform::parityBlocks(const int parity) const
   {
      if(parity == 0)
      {
         return this->mspSetup->evenBlocks();
      } else
      {
         return this->mspSetup->oddBlocks();
      }
   }

   inline fftw_plan SphereChebyshevFftwTransform::fPlan(const int inParity, const int outParity)
   {
      if(inParity == 0 && inParity == outParity)
      {
         return this->mFEPlan;
      } else if(inParity == 1 && inParity == outParity)
      {
         return this->mFBOPlan;
      } else if(inParity == 0)
      {
         return this->mFEOPlan;
      } else
      {
         return this->mFBOEPlan;
      }
   }

   inline fftw_plan SphereChebyshevFftwTransform::bPlan(const int inParity, const int outParity)
   {
      if(inParity == 0 && inParity == outParity)
      {
         return this->mBEPlan;
      } else if(inParity == 1 && inParity == outParity)
      {
         return this->mFBOPlan;
      } else if(inParity == 0)
      {
         return this->mBEOPlan;
      } else
      {
         return this->mFBOEPlan;
      }
   }

   inline SparseMatrix& SphereChebyshevFftwTransform::intgOp(const IntegratorType::Id integrator, const int parity)
   {
      if(parity == 0)
      {
         return this->mIntgOp.find(integrator)->second.first;
      } else
      {
         return this->mIntgOp.find(integrator)->second.second;
      }
   }

   inline SparseMatrix& SphereChebyshevFftwTransform::regOp(const RegularityType::Id reg, const int parity)
   {
      if(parity == 0)
      {
         return this->mRegOp.find(reg)->second.first;
      } else
      {
         return this->mRegOp.find(reg)->second.second;
      }
   }

   inline SparseMatrix& SphereChebyshevFftwTransform::solveOp(const ProjectorType::Id projector, const int parity)
   {
      if(parity == 0)
      {
         return this->mSolveOp.find(projector)->second.first;
      } else
      {
         return this->mSolveOp.find(projector)->second.second;
      }
   }

   inline SparseMatrix& SphereChebyshevFftwTransform::regSolveOp(const RegularityType::Id reg, const int parity)
   {
      if(parity == 0)
      {
         return this->mRegSolveOp.find(reg)->second.first;
      } else
      {
         return this->mRegSolveOp.find(reg)->second.second;
      }
   }

   inline Solver::SparseSelector<SparseMatrix>::Type& SphereChebyshevFftwTransform::solver(const ProjectorType::Id projector, const int parity)
   {
      if(parity == 0)
      {
         return *this->mSolver.find(projector)->second.first;
      } else
      {
         return *this->mSolver.find(projector)->second.second;
      }
   }

   inline Solver::SparseSelector<SparseMatrix>::Type& SphereChebyshevFftwTransform::regSolver(const RegularityType::Id reg, const int parity)
   {
      if(parity == 0)
      {
         return *this->mRegSolver.find(reg)->second.first;
      } else
      {
         return *this->mRegSolver.find(reg)->second.second;
      }
   }

   inline Matrix& SphereChebyshevFftwTransform::tmpInS(const int parity)
   {
      if(parity == 0)
      {
         return this->mTmpInSE;
      } else
      {
         return this->mTmpInSO;
      }
   }

   inline Matrix& SphereChebyshevFftwTransform::tmpOutS(const int parity)
   {
      if(parity == 0)
      {
         return this->mTmpOutSE;
      } else
      {
         return this->mTmpOutSO;
      }
   }

   inline int SphereChebyshevFftwTransform::regularitySize(const RegularityType::Id reg, const int parity) const
   {
      assert(this->mRegularitySize.count(reg) > 0);
   
      if(parity == 0)
      {
         return this->mRegularitySize.find(reg)->second.first;
      } else
      {
         return this->mRegularitySize.find(reg)->second.second;
      }
   }

   inline int SphereChebyshevFftwTransform::flipsParity(const ProjectorType::Id projector) const
   {
      assert(this->mProjectorFlips.count(projector) > 0);

      return this->mProjectorFlips.find(projector)->second;
   }

   inline int SphereChebyshevFftwTransform::flipsParity(const IntegratorType::Id integrator) const
   {
      assert(this->mIntegratorFlips.count(integrator) > 0);

      return this->mIntegratorFlips.find(integrator)->second;
   }

}
}

#endif // SPHERECHEBYSHEVFFTWTRANSFORM_HPP
