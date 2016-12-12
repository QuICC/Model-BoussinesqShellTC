/** 
 * @file BoussinesqFPlaneNHBGETransport.cpp
 * @brief Source of the implementation of the transport equation in the F-plane NHBGE model
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
#include "Equations/Asymptotics/FPlaneNHBGE/Boussinesq/BoussinesqFPlaneNHBGETransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqFPlaneNHBGETransport::BoussinesqFPlaneNHBGETransport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqFPlaneNHBGETransport::~BoussinesqFPlaneNHBGETransport()
   {
   }

   void BoussinesqFPlaneNHBGETransport::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 1, true, false);
   }

   void BoussinesqFPlaneNHBGETransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\theta\f$
      ///
      Physical::StreamAdvection<FieldComponents::Physical::X,FieldComponents::Physical::Y>::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);

      ///
      /// Computation of the mean temperature feedback
      ///
      rNLComp.addData((this->scalar(PhysicalNames::VELOCITYZ).dom(0).phys().data().array()*this->scalar(PhysicalNames::DZ_MEANTEMPERATURE).dom(0).phys().data().array()).matrix());
   }

   void BoussinesqFPlaneNHBGETransport::setRequirements()
   {
      // Set temperature as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, true));

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, true));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, false, true, false));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::DZ_MEANTEMPERATURE, FieldRequirement(true, false, true, false));

      // Gradient does not require Z component
      ArrayB   comps = ArrayB::Constant(3, true);
      comps(0) = false;
      std::map<FieldComponents::Spectral::Id,ArrayB>  gradComps;
      gradComps.insert(std::make_pair(FieldComponents::Spectral::SCALAR, comps));

      // Update temperature gradient requirements
      this->updateFieldRequirements(PhysicalNames::TEMPERATURE).updateGradient(gradComps);

      // Update streamfunction gradient requirements
      this->updateFieldRequirements(PhysicalNames::STREAMFUNCTION).updateGradient(gradComps);
   }

}
}
