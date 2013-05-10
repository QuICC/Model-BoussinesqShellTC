/** \file RayleighBenardTransport.cpp
 *  \brief Source of the implementation of the transport equation in the Rayleigh-Benard model
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
#include "Equations/Asymptotics/RayleighBenard/RayleighBenardTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RayleighBenardTransport::RayleighBenardTransport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   RayleighBenardTransport::~RayleighBenardTransport()
   {
   }

   void RayleighBenardTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T}\f$
      ///
   }

   void RayleighBenardTransport::setRequirements()
   {
      // Set temperature as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, true));
   }

   void RayleighBenardTransport::setCoupling()
   {
   }

   void RayleighBenardTransport::setSpectralMatrices(const IEvolutionEquation::BcEqMapType& bcIds, const std::map<PhysicalNames::Id, IEvolutionEquation::BcEqMapType>& cbcIds)
   {
   }
}
}
