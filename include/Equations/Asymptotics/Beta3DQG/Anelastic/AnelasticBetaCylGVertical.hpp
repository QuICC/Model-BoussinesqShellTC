/**
 * @file AnelasticBetaCylGVertical.hpp
 * @brief Implementation of the vertical velocity equation for the Anelastic beta model with cylindrical gravity 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ANELASTICBETACYLGVERTICAL_HPP
#define ANELASTICBETACYLGVERTICAL_HPP

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
    * @brief Implementation of the vertical velocity equation for the Anelastic beta model with cylindrical gravity
    */
   class AnelasticBetaCylGVertical: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         AnelasticBetaCylGVertical(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~AnelasticBetaCylGVertical();

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
         virtual DecoupledZSparse operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const;

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
         virtual void setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix &mat) const;

         /**
          * @brief Set the explicit linear matrix operator
          */
         virtual void setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const;

      private:
   };

   /**
    * @brief Get the quasi-inverse matrix operator
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    */
   void quasiInverseBlock(const AnelasticBetaCylGVertical& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat);

   /**
    * @brief Get the time matrix block
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param k       Wave number k
    */
   void timeBlock(const AnelasticBetaCylGVertical& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat k);

   /**
    * @brief Get the linear matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void linearBlock(const AnelasticBetaCylGVertical& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

   /**
    * @brief Get the boundary condition matrix block on given field
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @param k       Wave number k
    */
   void boundaryBlock(const AnelasticBetaCylGVertical& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k);

}
}

#endif // ANELASTICBETACYLGVERTICAL_HPP
