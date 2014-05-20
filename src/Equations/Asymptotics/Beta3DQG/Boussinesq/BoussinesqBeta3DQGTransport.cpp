/** 
 * @file BoussinesqBeta3DQGTransport.cpp
 * @brief Source of the implementation of the transport equation in the Beta 3DQG model
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBeta3DQGTransport::BoussinesqBeta3DQGTransport(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName,spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGTransport::~BoussinesqBeta3DQGTransport()
   {
   }

   void BoussinesqBeta3DQGTransport::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, true, false);
   }

   void BoussinesqBeta3DQGTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T}\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);
   }

   void BoussinesqBeta3DQGTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, true));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

}
}
