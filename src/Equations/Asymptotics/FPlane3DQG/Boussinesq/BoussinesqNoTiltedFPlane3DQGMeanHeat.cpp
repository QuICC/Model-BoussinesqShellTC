/** 
 * @file BoussinesqNoTiltedFPlane3DQGMeanHeat.cpp
 * @brief Source of the implementation of the mean heat computation in the tilted F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqNoTiltedFPlane3DQGMeanHeat.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqNoTiltedFPlane3DQGMeanHeat::BoussinesqNoTiltedFPlane3DQGMeanHeat(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqNoTiltedFPlane3DQGMeanHeat::~BoussinesqNoTiltedFPlane3DQGMeanHeat()
   {
   }

   void BoussinesqNoTiltedFPlane3DQGMeanHeat::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, true);
   }

   void BoussinesqNoTiltedFPlane3DQGMeanHeat::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get paramters
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat eta3 = std::cos((Math::PI/180.)*this->eqParams().nd(NonDimensional::THETA));
      MHDFloat eta2 = std::sin((Math::PI/180.)*this->eqParams().nd(NonDimensional::THETA));

      /// 
      /// Computation:
      ///   \f$ Pr w \theta\f$
      ///
      rNLComp.setData((eta3*Pr*this->scalar(PhysicalNames::NO_VELOCITYZ).dom(0).phys().data().array()*this->scalar(PhysicalNames::TEMPERATURE).dom(0).phys().data().array()).matrix());
      rNLComp.subData((eta2*Pr*this->scalar(PhysicalNames::NO_STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::X).data().array()*this->scalar(PhysicalNames::TEMPERATURE).dom(0).phys().data().array()).matrix());
   }

   Datatypes::SpectralScalarType::PointType BoussinesqNoTiltedFPlane3DQGMeanHeat::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iY) == 0)
      {
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iZ,iY) == 0)
         {
            if(iX == 0)
            {
               return Datatypes::SpectralScalarType::PointType(-1.0);
            }
         }
      } 

      return Datatypes::SpectralScalarType::PointType(0.0);
   }

   void BoussinesqNoTiltedFPlane3DQGMeanHeat::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::DZ_MEANTEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add mean temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::DZ_MEANTEMPERATURE, FieldRequirement(true, true, true, false));

      // Add temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, false, true, false));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::NO_STREAMFUNCTION, FieldRequirement(true, false, false, true));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::NO_VELOCITYZ, FieldRequirement(true, false, true, false));

      // Gradient does not require Z component
      ArrayB   comps = ArrayB::Constant(3, true);
      comps(0) = false;
      std::map<FieldComponents::Spectral::Id,ArrayB>  gradComps;
      gradComps.insert(std::make_pair(FieldComponents::Spectral::SCALAR, comps));

      // Update streamfunction gradient requirements
      this->updateFieldRequirements(PhysicalNames::NO_STREAMFUNCTION).updateGradient(gradComps);
   }

}
}
