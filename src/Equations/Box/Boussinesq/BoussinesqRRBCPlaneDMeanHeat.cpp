/** 
 * @file BoussinesqRRBCPlaneDMeanHeat.cpp
 * @brief Source of the implementation of the meant heat computation for rotating Rayleigh-Benard convection in a plane layer (toroidal/poloidal formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Box/Boussinesq/BoussinesqRRBCPlaneDMeanHeat.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/VelocityAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqRRBCPlaneDMeanHeat::BoussinesqRRBCPlaneDMeanHeat(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRRBCPlaneDMeanHeat::~BoussinesqRRBCPlaneDMeanHeat()
   {
   }

   void BoussinesqRRBCPlaneDMeanHeat::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, false, true);
   }

   Datatypes::SpectralScalarType::PointType BoussinesqRRBCPlaneDMeanHeat::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
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

   void BoussinesqRRBCPlaneDMeanHeat::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::DZ_MEANTEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::BEFORE);

      // Add mean temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::DZ_MEANTEMPERATURE, FieldRequirement(true, true, false, false));
   }

}
}
