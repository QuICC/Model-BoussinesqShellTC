/** \file RayleighBenardVertical.cpp
 *  \brief Source of the implementation of the vertical velocity equation in the Rayleigh-Benard model
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
#include "Equations/Asymptotics/RayleighBenard/RayleighBenardVertical.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RayleighBenardVertical::RayleighBenardVertical(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   RayleighBenardVertical::~RayleighBenardVertical()
   {
   }

   void RayleighBenardVertical::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ ??? \f$
      ///
   }

   void RayleighBenardVertical::setRequirements()
   {
      // Set vertical velocity as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));
   }

   void RayleighBenardVertical::setCoupling()
   {
   }
 
   void RayleighBenardVertical::setSpectralMatrices(const IEvolutionEquation::BcEqMapType& bcIds, const std::map<PhysicalNames::Id, IEvolutionEquation::BcEqMapType>& cbcIds)
   {
   }
}
}
