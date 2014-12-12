/**
 * @file CartesianExactVectorState.hpp
 * @brief Implementation of the equation to generate exact vector state in a cartesian geometry 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIANEXACTVECTORSTATE_HPP
#define CARTESIANEXACTVECTORSTATE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact vectors state in a cartesian geometry
    */
   class CartesianExactVectorState: public IVectorEquation
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
            // TTT states
            POLYPOLYPOLY = 10,  // Polynomial, Polynomial, Polynomial
            // TFT states
            POLYCOSPOLY = 20, // Polynomial, Cosine, Polynomial
            POLYSINPOLY,      // Polynomial, Sine, Polynomial
            // TFF states
            POLYCOSCOS = 30,  // Polynomial, Cosine, Cosine
            POLYSINSIN,       // Polynomial, Sine, Sine
            POLYSINCOS,       // Polynomial, Sine, Cosine
            POLYCOSSIN,       // Polynomial, Cosine, Sine
            // FFF states
            COSCOSCOS = 50,   // Cosine, Cosine, Cosine
            SINSINSIN,        // Sine, Sine, Sine
            COSCOSSIN,        // Cosine, Cosine, Sine
            SINSINCOS,        // Sine, Sine, Cosine
            COSSINSIN,        // Cosine, Sine, Sine
            SINCOSCOS,        // Sine, Cosine, Cosine
            COSSINCOS,        // Cosine, Sine, Cosine
            SINCOSSIN,        // Sine, Cosine, Sine
            // Special
            SPECIAL1 = 1000,  // Kind of a place holder for tests
            SPECIAL2,         // Kind of a place holder for tests
            SPECIAL3,         // Kind of a place holder for tests
         };

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         CartesianExactVectorState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~CartesianExactVectorState();

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
         void setStateType(const FieldComponents::Physical::Id compId, const CartesianExactVectorState::StateTypeId id);

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

      private:
         /**
          * @brief Compute even periodic mode
          */
         MHDFloat cos(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta) const;

         /**
          * @brief Compute odd periodic mode
          */
         MHDFloat sin(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta) const;

         /**
          * @brief Compute polynomial mode
          */
         MHDFloat poly(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x) const;

         /**
          * @brief Type of the state to generate
          */
         std::map<FieldComponents::Physical::Id,StateTypeId> mTypeId;

         /**
          * @brief Amplitude of the state
          */
         std::map<FieldComponents::Physical::Id,Array> mModeA;

         /**
          * @brief Mode number of the state (wave number of polynomial order)
          */
         std::map<FieldComponents::Physical::Id,Array> mModeK;
   };

   /// Typedef for a shared CartesianExactVectorState
   typedef SharedPtrMacro<CartesianExactVectorState> SharedCartesianExactVectorState;

}
}

#endif // CARTESIANEXACTVECTORSTATE_HPP
