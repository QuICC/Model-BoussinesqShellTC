/** 
 * @file BoussinesqFPlane3DQGMeanHeat.cpp
 * @brief Source of the implementation of the mean heat computation in the F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqFPlane3DQGMeanHeat.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "SpectralOperators/Tools/SpectralBoxTools.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

#include <iostream>
namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqFPlane3DQGMeanHeat::BoussinesqFPlane3DQGMeanHeat(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName,spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqFPlane3DQGMeanHeat::~BoussinesqFPlane3DQGMeanHeat()
   {
   }

   void BoussinesqFPlane3DQGMeanHeat::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::DIAGNOSTIC, 0, false, false, true);
   }

   void BoussinesqFPlane3DQGMeanHeat::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\nabla^2_{\perp}\psi\f$
      ///
      //Physical::StreamAdvection::set(rNLComp, this->unknown().dom(0).grad(), this->scalar(PhysicalNames::VORTICITYZ).dom(0).grad(), 1.0);
   }

   Datatypes::SpectralScalarType::PointType BoussinesqFPlane3DQGMeanHeat::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iY) == 0)
      {
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iZ,iY) == 0)
         {
            if(iX == 1)
            {
               return Datatypes::SpectralScalarType::PointType(-0.5);
            }
         }
      } 

      return Datatypes::SpectralScalarType::PointType(0.0);
   }

   void BoussinesqFPlane3DQGMeanHeat::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::MEANTEMPERATURE);

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::MEANTEMPERATURE, FieldRequirement(true, true, false, false));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, false));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, false));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, false, false));
   }

}
}
