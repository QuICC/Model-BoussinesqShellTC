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
#include <functional>

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

         /// Typedef for the function signature of an operator builder
         typedef std::function<void(AssociatedLegendreFlyTransform*, const int, const int, const int, const bool)> OperatorBuilder;

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
          * @brief Compute integration
          */
         void setIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, OperatorBuilder builder);

         /**
          * @brief Compute integration with f(l) prefactor
          */
         void setMultLIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const Array& mult, OperatorBuilder builder);

         /**
          * @brief Compute projection
          */
         void setProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, OperatorBuilder builder);

         /**
          * @brief Compute projection with f(l) prefactor
          */
         void setMultLProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, const Array& mult, OperatorBuilder builder);

         /**
          * @brief Build Plm operator for m+start to m+start+nL-1 
          */
         void buildPlm(const int idx, const int start, const int nL, const bool isWeighted);

         /**
          * @brief Build \partial_\theta Plm operator for m+start to m+start+nL-1 
          */
         void buildDPlm(const int idx, const int start, const int nL, const bool isWeighted);

         /**
          * @brief Build \frac{1}{\sin(\theta)} Plm operator for m+start to m+start+nL-1 
          */
         void buildSin_1Plm(const int idx, const int start, const int nL, const bool isWeighted);

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
          * @brief Storage for operator
          */
         Matrix mOp;

         /**
          * @brief Temporary storage for running recurrence A 
          */
         Eigen::Matrix<MHDFloat, Eigen::Dynamic, 2> mTmpA;

         /**
          * @brief Temporary storage for running recurrence B 
          */
         Eigen::Matrix<MHDFloat, Eigen::Dynamic, 2> mTmpB;
   };

}
}

#endif // ASSOCIATEDLEGENDREFLYTRANSFORM_HPP
