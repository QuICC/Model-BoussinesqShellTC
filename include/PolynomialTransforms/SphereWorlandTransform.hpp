/** 
 * @file SphereWorlandTransform.hpp
 * @brief Implementation of the Worland transform in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHEREWORLANDTRANSFORM_HPP
#define SPHEREWORLANDTRANSFORM_HPP

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

namespace QuICC {

namespace Transform {

   /**
    * @brief Simple struct holding details about Worland transform
    */
   struct SphereWorlandIds {

      /**
       * @brief Simple struct holding the projector IDs
       */
      struct Projectors
      {
         /** Enum of projector IDs
          *    - PROJ: projection
          *    - DIVR: 1/r
          *    - DIFF: D
          *    - DIFFR: D r
          *    - DIVRDIFFR: 1/r D r
          *    - SLAPL: spherical laplacian: D^2 + 2/r D - l(l+1)/r^2
          *    - ENERGY_PROJ:
          *    - ENERGY_DIFFR:
          */
         enum Id {PROJ, DIVR, DIFF, DIFFR, DIVRDIFFR, SLAPL, ENERGY_PROJ, ENERGY_DIFFR};
      };

      /**
       * @brief Simple struct holding the integrator IDs
       */
      struct Integrators
      {
         /**
          * Enum of integrator IDs
          *    - INTG: integration
          *    - INTGR: integration of r
          *    - INTGQ4: integration of QST Q component for Poloidal NL (4th order equation)
          *    - INTGS4: integration of QST S component for Poloidal NL (4th order equation)
          *    - INTGT: integration of QST T component for Toroidal NL (2nd order equation)
          *    - INTGQ2: integration of QST Q component for Poloidal NL (2nd order equation)
          *    - INTGS2: integration of QST S component for Poloidal NL (2nd order equation)
          *    - ENERGY_INTG: compute definite energy integral
          *    - ENERGY_R2: compute definite energy integral with spherical radial weight
          */
         enum Id {INTG, INTGR, INTGQ4, INTGS4, INTGT, INTGQ2, INTGS2, ENERGY_INTG, ENERGY_R2};
      };

   };

   /**
    * @brief Implementation of the Worland transform in a sphere
    */ 
   class SphereWorlandTransform
   {
      public:
         /// Typedef for the configuration class
         typedef PolySetup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedPolySetup SharedSetupType;

         /// Typedef for the Projector type
         typedef SphereWorlandIds::Projectors ProjectorType;

         /// Typedef for the Integrator type
         typedef SphereWorlandIds::Integrators IntegratorType;

         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size); 

         /**
          * @brief Constructor
          */
         SphereWorlandTransform();

         /**
          * @brief Destructor
          */
         ~SphereWorlandTransform();

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

         /**
          * @brief Compute quadrature integration for energy integral
          *
          * Compute quadrature integration from energy integral
          *
          * @param spectrum   Output energy spectrum
          * @param specVal    Input physical values
          * @param integrator Integrator to use
          */
         void integrate_energy(Array& spectrum, const MatrixZ& specVal, ProjectorType::Id projector, IntegratorType::Id integrator);

     #ifdef QUICC_STORAGEPROFILE
         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
     #endif // QUICC_STORAGEPROFILE
         
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
          * @brief Compute energy integration with vector of operators
          */
         void setEnergyIntegrator(Array& spectrum, const MatrixZ& specVal, const std::vector<Matrix>& projOps, const std::vector<Matrix>& intgOps);

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
   };

}
}

#endif // SPHEREWORLANDTRANSFORM_HPP
