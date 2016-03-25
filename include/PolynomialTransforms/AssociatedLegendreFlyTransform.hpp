/** 
 * @file AssociatedLegendreFlyTransform.hpp
 * @brief Implementation of the associated Legendre on-the-fly transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ASSOCIATEDLEGENDREFLYTRANSFORM_HPP
#define ASSOCIATEDLEGENDREFLYTRANSFORM_HPP

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
#include "Enums/Arithmetics.hpp"
#include "Enums/NonDimensional.hpp"
#include "PolynomialTransforms/PolySetup.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Simple struct holding details about associated Legendre transform
    */
   struct AssociatedLegendreIds {

      /**
       * @brief Simple struct holding the projector IDs
       */
      struct Projectors
      {
         /// Enum of projector IDs
         enum Id {PROJ, PROJLL, DIFF, DIFFLL, DIVSIN, DIVSINLL, DIVSINDIFFSIN};
      };

      /**
       * @brief Simple struct holding the integrator IDs
       */
      struct Integrators
      {
         /// Enum of integrator IDs
         enum Id {INTG, INTGDIVLL, INTGLL, INTGLL2, INTGDIVSIN, INTGDIVLLDIVSIN, INTGLLDIVSIN, INTGDIFF, INTGDIVLLDIFF, INTGLLDIFF};
      };

   };

   /**
    * @brief Implementation of the associated Legendre on-the-fly transform
    */ 
   class AssociatedLegendreFlyTransform
   {
      public:
         /// Typedef for the configuration class
         typedef PolySetup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedPolySetup SharedSetupType;

         /// Typedef for the Projector type
         typedef AssociatedLegendreIds::Projectors ProjectorType;

         /// Typedef for the Integrator type
         typedef AssociatedLegendreIds::Integrators IntegratorType;

         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size); 

         /**
          * @brief Constructor
          */
         AssociatedLegendreFlyTransform();

         /**
          * @brief Destructor
          */
         ~AssociatedLegendreFlyTransform();

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
          * @param arithId    Arithmetic operation to perform
          */
         void integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, IntegratorType::Id integrator, Arithmetics::Id arithId);

         /**
          * @brief Compute polynomial projection
          *
          * @param rPhysVal   Output physical values
          * @param specVal    Input spectral coefficients
          * @param projector  Projector to use
          * @param arithId    Arithmetic operation to perform
          */
         void project(MatrixZ& rPhysVal, const MatrixZ& specVal, ProjectorType::Id projector, Arithmetics::Id arithId);

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
         void setIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal);

         /**
          * @brief Compute projection with vector of operators
          */
         void setProjector(MatrixZ& rPhysVal, const MatrixZ& specVal);

         /**
          * @brief Build Plm operator for m+start to m+start+nL-1 
          */
         void buildPlm(const int idx, const int start, const int nL);

         /**
          * @brief Build weighted Plm operator for m+start to m+start+nL-1 
          */
         void buildWeightedPlm(const int idx, const int start, const int nL);

         /**
          * @brief Max size of the operator
          */
         const int mcMaxOpCols;

         /**
          * @brief Size of the operator
          */
         int mOpCols;

         /**
          * @brief Storage for the quadrature points x = [-1, 1]
          */
         Array mXGrid;

         /**
          * @brief Storage for the quadrature points th = [0, pi]
          */
         Array mThGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         Array mWeights;

         /**
          * @brief Storage for the l(l+1) factor
          */
         Array mLl;

         /**
          * @brief Storage for the (l(l+1))^2 factor
          */
         Array mLl2;

         /**
          * @brief Storage for the 1/l(l+1) factor
          */
         Array mDivLl;

         /**
          * @brief Polynomial setup object providing the sizes
          */
         SharedPolySetup    mspSetup;

         /**
          * @brief Storage for A recurrence term
          */
         Matrix mOp;
   };

}
}

#endif // ASSOCIATEDLEGENDREFLYTRANSFORM_HPP
