/** 
 * @file CylinderWorlandTransform.hpp
 * @brief Implementation of the Worland transform in a cylinder
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CYLINDERWORLANDTRANSFORM_HPP
#define CYLINDERWORLANDTRANSFORM_HPP

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

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"
#include "Enums/NonDimensional.hpp"
#include "PolynomialTransforms/PolySetup.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "SparseSolvers/SparseLinearSolverTools.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Simple struct holding details about Worland transform
    */
   struct CylinderWorlandIds {

      /**
       * @brief Simple struct holding the projector IDs
       */
      struct Projectors
      {
         /** Enum of projector IDs
          *    - PROJ: projection
          *    - DIVRM0: 1/r, 0 for m = 0
          *    - DIFF: D
          *    - DIVRDIFFR: 1/r D r
          *    - LAPLH: horizontal cylindrical laplacian: D^2 + 1/r D - m^2/r^2
          *    - DIVRLAPLHM0: 1/r horizontal cylindrical laplacian: D^2 + 1/r D - m^2/r^2, 0 for m = 0
          *    - DIFFLAPLH: radial derivative of horizontal cylindrical laplacian: D(D^2 + 1/r D - m^2/r^2)
          */
         enum Id {PROJ, DIVRM0, DIFF, DIVRDIFFR, LAPLH, DIVRLAPLHM0, DIFFLAPLH, INVLAPLH};
      };

      /**
       * @brief Simple struct holding the integrator IDs
       */
      struct Integrators
      {
         /**
          * Enum of integrator IDs
          *    - INTG: integration
          *    - INTGI4DIVRM0: integration, zero for m = 0
          *    - INTGI4DIVRDIFFR: integration
          *    - INTGI6DIVR: integration, 0 for m = 0
          */
         enum Id {INTG, INTGI4DIVRM0, INTGI4DIVRDIFFR, INTGI6DIVRM0, INTGI6DIVRDIFFR, INTGI6LAPLH, INTGINVLAPLH};
      };

   };

   /**
    * @brief Implementation of the Worland transform in a sphere
    */ 
   class CylinderWorlandTransform
   {
      public:
         /// Typedef for the configuration class
         typedef PolySetup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedPolySetup SharedSetupType;

         /// Typedef for the Projector type
         typedef CylinderWorlandIds::Projectors ProjectorType;

         /// Typedef for the Integrator type
         typedef CylinderWorlandIds::Integrators IntegratorType;

         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size); 

         /**
          * @brief Constructor
          */
         CylinderWorlandTransform();

         /**
          * @brief Destructor
          */
         ~CylinderWorlandTransform();

         /**
          * @brief Initialise the polynomial transform (matrices, weights, grid, etc)
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
         const Array& meshGrid() const; 

         /**
          * @brief Compute quadrature integration
          *
          * @param rSpecVal   Output spectral coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          */
         void integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute polynomial projection
          *
          * @param rPhysVal   Output physical values
          * @param specVal    Input spectral coefficients
          * @param projector  Projector to use
          */
         void project(MatrixZ& rPhysVal, const MatrixZ& specVal, ProjectorType::Id projector);

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief Initialise the operators
          */
         void initOperators();

         /**
          * @brief Compute integration with vector of operators
          */
         void setIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const std::vector<Matrix>& ops);

         /**
          * @brief Compute projection with vector of operators
          */
         void setProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, const std::vector<Matrix>& ops);

         /**
          * @brief Storage for the quadrature points x = [-1, 1]
          */
         Array mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         Array mWeights;

         /**
          * @brief Polynomial setup object providing the sizes
          */
         SharedPolySetup    mspSetup;

         /**
          * @brief Storage for the projector operators 
          */
         std::map<ProjectorType::Id,std::vector<Matrix> >  mProjOp;

         /**
          * @brief Storage for the integrator operators 
          */
         std::map<IntegratorType::Id,std::vector<Matrix> >  mIntgOp;

         /**
          * @brief Storage for the sparse solver matrices
          */
         std::map<ProjectorType::Id, std::vector<SparseMatrix> > mSolveOp;

         /**
          * @brief Storage for the sparse triangular solvers
          */
         std::map<ProjectorType::Id, std::vector<SharedPtrMacro<Solver::SparseTriSelector<SparseMatrix>::Type> > > mTriSolver;
   };

}
}

#endif // CYLINDERWORLANDTRANSFORM_HPP
