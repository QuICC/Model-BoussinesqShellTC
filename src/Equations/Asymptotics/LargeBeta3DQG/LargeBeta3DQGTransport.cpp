/** \file LargeBeta3DQGTransport.cpp
 *  \brief Source of the implementation of the transport equation in the large 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Asymptotics/LargeBeta3DQG/LargeBeta3DQGTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   LargeBeta3DQGTransport::LargeBeta3DQGTransport(SharedIEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   LargeBeta3DQGTransport::~LargeBeta3DQGTransport()
   {
   }

   void LargeBeta3DQGTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T}\f$
      ///
   }

   void LargeBeta3DQGTransport::setRequirements()
   {
      // Set temperature as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, true));
   }

   void LargeBeta3DQGTransport::setCoupling()
   {
   }

   void LargeBeta3DQGTransport::setSpectralMatrices(const IEvolutionEquation::BcEqMapType& bcIds, const std::map<PhysicalNames::Id, IEvolutionEquation::BcEqMapType>& cbcIds)
   {
   }
}
}