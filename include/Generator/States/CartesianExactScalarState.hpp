/**
 * @file CartesianExactScalarState.hpp
 * @brief Implementation of the equation to generate exact scalar states in a cartesian geometries 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIANEXACTSCALARSTATE_HPP
#define CARTESIANEXACTSCALARSTATE_HPP

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
    * @brief Implementation of the equation to generate exact scalar state in cartesian geometries
    */
   class CartesianExactScalarState: public IScalarEquation
   {
      public:
         /**
          * @brief Enums for the avaialable exact states
          */
         enum StateTypeId {
            CONSTANT = 0,  // All constant
            POLYPOLYPOLY = 10,  // Polynomial, Polynomial, Polynomial
            POLYCOSPOLY = 20, // Polynomial, Cosine, Polynomial
            POLYSINPOLY,      // Polynomial, Sine, Polynomial
            POLYCOSCOS = 30,  // Polynomial, Cosine, Cosine
            POLYSINSIN,       // Polynomial, Sine, Sine
            POLYSINCOS,       // Polynomial, Sine, Cosine
            POLYCOSSIN,       // Polynomial, Cosine, Sine
            COSCOSCOS = 40,   // Cosine, Cosine, Cosine
            SINSINSIN,        // Sine, Sine, Sine
            COSCOSSIN,        // Cosine, Cosine, Sine
            SINSINCOS,        // Sine, Sine, Cosine
            COSSINSIN,        // Cosine, Sine, Sine
            SINCOSCOS,        // Sine, Cosine, Cosine
            COSSINCOS,        // Cosine, Sine, Cosine
            SINCOSSIN,        // Sine, Cosine, Sine
         };

         /**
          * @brief Simple constructor
          *
          * @param pyName     Python script name
          * @param spEqParams Shared equation parameters
          */
         CartesianExactScalarState(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~CartesianExactScalarState();

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
         void setStateType(const CartesianExactScalarState::StateTypeId id);

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

   /// Typedef for a shared CartesianExactScalarState
   typedef SharedPtrMacro<CartesianExactScalarState> SharedCartesianExactScalarState;

}
}

#endif // CARTESIANEXACTSCALARSTATE_HPP
