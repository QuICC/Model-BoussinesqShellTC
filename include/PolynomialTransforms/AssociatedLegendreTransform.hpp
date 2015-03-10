/** 
 * @file AssociatedLegendreTransform.hpp
 * @brief Implementation of the associated Legendre transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ASSOCIATEDLEGENDRETRANSFORM_HPP
#define ASSOCIATEDLEGENDRETRANSFORM_HPP

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
    * @brief Simple struct holding details about ChebyshevFFT transform
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
    * @brief Implementation of the FFTW transform for a Chebyshev expansion
    */ 
   class AssociatedLegendreTransform
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
         AssociatedLegendreTransform();

         /**
          * @brief Destructor
          */
         ~AssociatedLegendreTransform();

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
         void setIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const std::vector<Matrix>& ops);

         /**
          * @brief l(l+1) Compute integration with vector of operators
          */
         void setLlIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const std::vector<Matrix>& ops);

         /**
          * @brief Compute l^2(l+1)^2 integration with vector of operators
          */
         void setLl2Integrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const std::vector<Matrix>& ops);

         /**
          * @brief Compute 1/l(l+1) integration with vector of operators
          */
         void setDivLlIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const std::vector<Matrix>& ops);

         /**
          * @brief Compute projection with vector of operators
          */
         void setProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, const std::vector<Matrix>& ops);

         /**
          * @brief Compute  l(l+1) projection with vector of operators
          */
         void setLlProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, const std::vector<Matrix>& ops);

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
          * @brief Projector matrix
          */
         std::vector<Matrix>  mProj;

         /**
          * @brief Derivative matrix
          */
         std::vector<Matrix>  mDiff;

         /**
          * @brief \f$1/\sin\theta\f$ projector matrix
          */
         std::vector<Matrix>  mDivSin;
   };

}
}

#endif // ASSOCIATEDLEGENDRETRANSFORM_HPP
