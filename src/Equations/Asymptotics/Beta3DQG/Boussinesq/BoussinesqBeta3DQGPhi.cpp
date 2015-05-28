/** 
 * @file BoussinesqBeta3DQGPhi.cpp
 * @brief Source of the implementation of the phi equation (D_z w) in the Beta 3DQG model
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPhi.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBeta3DQGPhi::BoussinesqBeta3DQGPhi(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGPhi::~BoussinesqBeta3DQGPhi()
   {
   }

   void BoussinesqBeta3DQGPhi::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqBeta3DQGPhi::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///
      Physical::StreamAdvection<>::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->scalar(PhysicalNames::VELOCITYZ).dom(0).grad(), 1.0);
   }

   void BoussinesqBeta3DQGPhi::setRequirements()
   {
      // Set vertical velocity as equation unknown
      this->setName(PhysicalNames::PHI);

      // Add phi requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::PHI, FieldRequirement(true, true, false, false));

      // Add D_z phi requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, false, true));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, true));
   }

}
}
