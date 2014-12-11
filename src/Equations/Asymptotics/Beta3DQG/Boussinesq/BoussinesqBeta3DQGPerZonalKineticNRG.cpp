/** 
 * @file BoussinesqBeta3DQGPerZonalKineticNRG.cpp
 * @brief Source of the implementation of the zonal kinetic energy caluclation in the periodic Beta 3DQG model
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerZonalKineticNRG.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBeta3DQGPerZonalKineticNRG::BoussinesqBeta3DQGPerZonalKineticNRG(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGPerZonalKineticNRG::~BoussinesqBeta3DQGPerZonalKineticNRG()
   {
   }

   void BoussinesqBeta3DQGPerZonalKineticNRG::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, true, false);
   }

   void BoussinesqBeta3DQGPerZonalKineticNRG::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation:
      ///   \f$ u \cdot u\f$
      ///
      rNLComp.setData((this->scalar(PhysicalNames::MEAN_VELOCITYY).dom(0).phys().data().array().pow(2)).matrix());
      rNLComp.addData((this->scalar(PhysicalNames::MEAN_VELOCITYZ).dom(0).phys().data().array().pow(2)).matrix());
   }

   void BoussinesqBeta3DQGPerZonalKineticNRG::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::ZONAL_KINETIC_ENERGY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add kinetic energy to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::ZONAL_KINETIC_ENERGY, FieldRequirement(true, true, false, false));

      // Add mean Y velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::MEAN_VELOCITYY, FieldRequirement(true, true, true, false));

      // Add mean Z velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::MEAN_VELOCITYZ, FieldRequirement(true, true, true, false));
   }

}
}
