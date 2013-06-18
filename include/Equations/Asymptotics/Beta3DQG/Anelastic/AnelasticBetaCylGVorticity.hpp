/** \file AnelasticBetaCylGVorticity.hpp
 *  \brief Implementation of the voriticity computation for the Anelastic beta model with cylindrical gravity
 */

#ifndef ANELASTICBETACYLGVORTICITY_HPP
#define ANELASTICBETACYLGVORTICITY_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the voriticity computation for the Anelastic beta model with cylindrical gravity
    */
   class AnelasticBetaCylGVorticity: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         AnelasticBetaCylGVorticity(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~AnelasticBetaCylGVorticity();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Generic operator row dispatcher
          */
         virtual DecoupledZSparse operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id comp, const int matIdx) const;

         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds);
         
      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the quasi inverse matrix operator
          */
         virtual void setQuasiInverse(SparseMatrix &mat) const;

         /**
          * @brief Set the explicit linear matrix operator
          */
         virtual void setExplicitLinearBlock(DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const;

      private:
   };

   /**
    * @brief Get the linear matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void linearBlock(const AnelasticBetaCylGVorticity& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

}
}

#endif // ANELASTICBETACYLGVORTICITY_HPP
