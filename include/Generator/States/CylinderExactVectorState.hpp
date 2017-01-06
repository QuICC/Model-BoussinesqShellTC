/**
 * @file CylinderExactVectorState.hpp
 * @brief Implementation of the equation to generate exact vector states in a cylinder 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CYLINDEREXACTVECTORSTATE_HPP
#define CYLINDEREXACTVECTORSTATE_HPP

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
#include "Equations/IVectorEquation.hpp"
#include "Generator/States/CylinderExactStateIds.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact vector state in a cylinder
    */
   class CylinderExactVectorState: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         CylinderExactVectorState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~CylinderExactVectorState();

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
         void setPhysicalType(const FieldComponents::Physical::Id compId, const CylinderExactStateIds::Id id);

         /**
          * @brief Set the state type id
          */
         void setSpectralType(const FieldComponents::Spectral::Id compId, const CylinderExactStateIds::Id id);

         /**
          * @brief Set the options for the solution states
          *
          * @param a1   Amplitude of the first direction
          * @param k1   Wave number of the first direction
          * @param a2   Amplitude of the second direction
          * @param k2   Wave number of the second direction
          * @param a3   Amplitude of the second direction
          * @param k3   Wave number of the second direction
          */
         void setModeOptions(const FieldComponents::Physical::Id compId, const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the nonliner integration components
          */
         virtual void setNLComponents();

      private:
         /**
          * @brief Type of the state to generate
          */
         std::map<FieldComponents::Physical::Id,CylinderExactStateIds::Id> mTypeId;

         /**
          * @brief Type of the spectral state to generate
          */
         std::map<FieldComponents::Spectral::Id,CylinderExactStateIds::Id> mSpecTypeId;

         /**
          * @brief Amplitude of the state
          */
         std::map<FieldComponents::Physical::Id,Array> mModeA;

         /**
          * @brief Mode number of the state (wave number of polynomial order)
          */
         std::map<FieldComponents::Physical::Id,Array> mModeK;
   };

   /// Typedef for a shared CylinderExactVectorState
   typedef SharedPtrMacro<CylinderExactVectorState> SharedCylinderExactVectorState;

}
}

#endif // CYLINDEREXACTVECTORSTATE_HPP
