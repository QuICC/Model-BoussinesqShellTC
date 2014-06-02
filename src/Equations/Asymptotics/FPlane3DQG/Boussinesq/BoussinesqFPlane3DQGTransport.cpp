/** 
 * @file BoussinesqFPlane3DQGTransport.cpp
 * @brief Source of the implementation of the transport equation in the F-plane 3DQG model
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
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqFPlane3DQGTransport::BoussinesqFPlane3DQGTransport(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName,spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqFPlane3DQGTransport::~BoussinesqFPlane3DQGTransport()
   {
   }

   void BoussinesqFPlane3DQGTransport::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 1, true, false, false);
   }

   void BoussinesqFPlane3DQGTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get parameters
      MHDFloat eta2 = std::sin((Math::PI/180.)*this->eqParams().nd(NonDimensional::THETA));
      MHDFloat eta3 = std::cos((Math::PI/180.)*this->eqParams().nd(NonDimensional::THETA));

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\thetaT\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0/eta3, FieldComponents::Physical::TWO, FieldComponents::Physical::THREE);

      ///
      /// Computation of the mean temperature feedback
      ///
      rNLComp.addData((2.0*(-eta3*this->scalar(PhysicalNames::VELOCITYZ).dom(0).phys().data().array() + eta2*this->unknown().dom(0).grad().comp(FieldComponents::Physical::TWO).data().array())*this->scalar(PhysicalNames::MEANTEMPERATURE).dom(0).grad().comp(FieldComponents::Physical::ONE).data().array()).matrix());
   }

   void BoussinesqFPlane3DQGTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, true));

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, true));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, false, true, false));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::MEANTEMPERATURE, FieldRequirement(true, false, false, true));
   }

}
}
