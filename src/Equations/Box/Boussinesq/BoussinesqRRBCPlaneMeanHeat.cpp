/** 
 * @file BoussinesqRRBCPlaneMeanHeat.cpp
 * @brief Source of the implementation of the meant heat computation for rotating Rayleigh-Benard convection in a plane layer (toroidal/poloidal formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Box/Boussinesq/BoussinesqRRBCPlaneMeanHeat.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRRBCPlaneMeanHeat::BoussinesqRRBCPlaneMeanHeat(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRBCPlaneMeanHeat::~BoussinesqRRBCPlaneMeanHeat()
   {
   }

   void BoussinesqRRBCPlaneMeanHeat::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::DIAGNOSTIC, 0, false, false);
   }

   void BoussinesqRRBCPlaneMeanHeat::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
   }

   void BoussinesqRRBCPlaneMeanHeat::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::MEAN_TEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::BEFORE);

      // Add mean temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::MEAN_TEMPERATURE, FieldRequirement(true, true, false, false));
   }

}
}
