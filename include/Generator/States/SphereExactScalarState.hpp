/**
 * @file SphereExactScalarState.hpp
 * @brief Implementation of the equation to generate exact scalar states in a sphere 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHEREEXACTSCALARSTATE_HPP
#define SPHEREEXACTSCALARSTATE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <vector>
#include <tuple>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Generator/States/SphereExactStateIds.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact scalar state in a sphere
    */
   class SphereExactScalarState: public IScalarEquation
   {
      public:
         /// Typedef to simplify notations for harmonic mode
         typedef std::map<std::pair<int, int>, std::map<int,MHDComplex> > HarmonicModeType;

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         SphereExactScalarState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~SphereExactScalarState();

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
         void setStateType(const SphereExactStateIds::Id id);

         /**
          * @brief Set the spectral state type id
          */
         void setSpectralType(const SphereExactStateIds::Id id);

         /**
          * @brief Set options for the harmonics states
          *
          * @param modes   List of harmonics with amplitude to create
          */
         void setHarmonicOptions(const HarmonicModeType& modes);

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
         SphereExactStateIds::Id mTypeId;

         /**
          * @brief Type of the spectral state to generate
          */
         SphereExactStateIds::Id mSpecTypeId;

         /**
          * @brief Storage for the list of spherical harmonic modes to generate
          */
         HarmonicModeType mSHModes;
   };

   /// Typedef for a shared SphereExactScalarState
   typedef SharedPtrMacro<SphereExactScalarState> SharedSphereExactScalarState;

}
}

#endif // SPHEREEXACTSCALARSTATE_HPP
