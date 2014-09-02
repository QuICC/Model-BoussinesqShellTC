/**
 * @file CylinderExactScalarState.hpp
 * @brief Implementation of the equation to generate exact scalar states in a cylinder 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CYLINDEREXACTSCALARSTATE_HPP
#define CYLINDEREXACTSCALARSTATE_HPP

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
    * @brief Implementation of the equation to generate exact scalar state in a cylinder
    */
   class CylinderExactScalarState: public IScalarEquation
   {
      public:
         /// Polynomial approximation to Cosine
         static const MHDFloat PCOS = 99999;

         /// Polynomial approximation to Sine
         static const MHDFloat PSIN = -99999;

         /**
          * @brief Enums for the avaialable exact states
          */
         enum StateTypeId {
            // Special states
            CONSTANT = 0,  // All constant
            // CFT states
            POLYCOSPOLY = 20, // Polynomial, Cosine, Polynomial
            POLYSINPOLY,      // Polynomial, Sine, Polynomial
         };

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         CylinderExactScalarState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~CylinderExactScalarState();

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
         void setStateType(const CylinderExactScalarState::StateTypeId id);

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
         void setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3);

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
          * @brief Compute even periodic mode
          */
         MHDFloat cos(const int idx, const MHDFloat theta) const;

         /**
          * @brief Compute odd periodic mode
          */
         MHDFloat sin(const int idx, const MHDFloat theta) const;

         /**
          * @brief Compute polynomial mode
          */
         MHDFloat poly(const int idx, const MHDFloat x) const;

         /**
          * @brief Type of the state to generate
          */
         StateTypeId mTypeId;

         /**
          * @brief Amplitude of the state
          */
         Array mModeA;

         /**
          * @brief Mode number of the state (wave number of polynomial order)
          */
         Array mModeK;
   };

   /// Typedef for a shared CylinderExactScalarState
   typedef SharedPtrMacro<CylinderExactScalarState> SharedCylinderExactScalarState;

}
}

#endif // CYLINDEREXACTSCALARSTATE_HPP
