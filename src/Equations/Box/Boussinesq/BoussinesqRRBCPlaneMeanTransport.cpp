/** 
 * @file BoussinesqRRBCPlaneMeanTransport.cpp
 * @brief Source of the implementation of the transport equation for rotating Rayleigh-Benard convection in a plane layer (toroidal/poloidal formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Box/Boussinesq/BoussinesqRRBCPlaneMeanTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRRBCPlaneMeanTransport::BoussinesqRRBCPlaneMeanTransport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRBCPlaneMeanTransport::~BoussinesqRRBCPlaneMeanTransport()
   {
   }

   void BoussinesqRRBCPlaneMeanTransport::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqRRBCPlaneMeanTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)\theta\f$
      ///
      Physical::VelocityAdvection<FieldComponents::Physical::X,FieldComponents::Physical::Y,FieldComponents::Physical::Z>::set(rNLComp, this->vector(PhysicalNames::VELOCITY).dom(0).phys(), this->unknown().dom(0).grad(), 1.0);

      ///
      /// Computation of the mean temperature feedback
      ///
      rNLComp.addData((this->vector(PhysicalNames::VELOCITY).dom(0).phys().comp(FieldComponents::Physical::Z).data().array()*this->scalar(PhysicalNames::MEAN_TEMPERATURE).dom(0).grad().comp(FieldComponents::Physical::Z).data().array()).matrix());
   }

   void BoussinesqRRBCPlaneMeanTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, true));

      // Add mean temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::MEAN_TEMPERATURE, FieldRequirement(true, true, false, true));

      // Add X velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false));

      // Mean temperature gradient only requires Z component
      ArrayB   comps = ArrayB::Constant(3, false);
      comps(0) = true;
      std::map<FieldComponents::Spectral::Id,ArrayB>  gradComps;
      gradComps.insert(std::make_pair(FieldComponents::Spectral::SCALAR, comps));

      // Update mean temperature gradient requirements
      this->updateFieldRequirements(PhysicalNames::MEAN_TEMPERATURE).updateGradient(gradComps);
   }

}
}
