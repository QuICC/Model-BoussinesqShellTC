/** \file RayleighBenardStreamfunction.cpp
 *  \brief Source of the implementation of the streamfunction equation in the Rayleigh-Benard model
 */

// Configuration includes
//

// System includes
//

// External includes
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Asymptotics/RayleighBenard/RayleighBenardStreamfunction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RayleighBenardStreamfunction::RayleighBenardStreamfunction(SharedEquationParameters spEqParams)
      : IScalarPEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   RayleighBenardStreamfunction::~RayleighBenardStreamfunction()
   {
   }

   void RayleighBenardStreamfunction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ ??? \f$
      ///
   }

   void RayleighBenardStreamfunction::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::STREAMFUNCTION);

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, true, true));
   }

   void RayleighBenardStreamfunction::setCoupling()
   {
   }

   void RayleighBenardStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarPEquation::timestepOutput(id, storage, matIdx, start);
   }

   void RayleighBenardStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarPEquation::timestepOutput(id, storage, matIdx, start);
   }
 
   void RayleighBenardStreamfunction::setSpectralMatrices(const IEvolutionEquation::BcEqMapType& bcIds, const std::map<PhysicalNames::Id, IEvolutionEquation::BcEqMapType>& cbcIds)
   {
   }
}
}
