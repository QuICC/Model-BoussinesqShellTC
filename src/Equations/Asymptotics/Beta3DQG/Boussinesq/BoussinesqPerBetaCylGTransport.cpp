/** \file BoussinesqPerBetaCylGTransport.cpp
 *  \brief Source of the implementation of the transport equation in the 3DQG beta model with periodic radius
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGSystem.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqPerBetaCylGTransport::BoussinesqPerBetaCylGTransport(SharedEquationParameters spEqParams)
      : IBoussinesqPerBetaCylGScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqPerBetaCylGTransport::~BoussinesqPerBetaCylGTransport()
   {
   }

   void BoussinesqPerBetaCylGTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T}\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);
   }

   void BoussinesqPerBetaCylGTransport::setRequirements()
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
