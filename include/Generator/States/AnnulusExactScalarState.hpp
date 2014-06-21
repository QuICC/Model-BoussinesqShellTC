/**
 * @file AnnulusExactScalarState.hpp
 * @brief Implementation of the equation to generate exact scalar states in an annulus 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ANNULUSEXACTSCALARSTATE_HPP
#define ANNULUSEXACTSCALARSTATE_HPP

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
    * @brief Implementation of the equation to generate exact scalar state in an annulus
    */
   class AnnulusExactScalarState: public IScalarEquation
   {
      public:
         /**
          * @brief Enums for the avaialable exact states
          */
         enum StateTypeId {
            CONSTANT,
         };

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         AnnulusExactScalarState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~AnnulusExactScalarState();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const;

         /**
          * @brief Compute the source term
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual Datatypes::SpectralScalarType::PointType sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set the state type id
          */
         void setStateType(const AnnulusExactScalarState::StateTypeId id);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
         /**
          * @brief Type of the state to generate
          */
         StateTypeId mTypeId;
   };

   /// Typedef for a shared AnnulusExactScalarState
   typedef SharedPtrMacro<AnnulusExactScalarState> SharedAnnulusExactScalarState;

}
}

#endif // ANNULUSEXACTSCALARSTATE_HPP
