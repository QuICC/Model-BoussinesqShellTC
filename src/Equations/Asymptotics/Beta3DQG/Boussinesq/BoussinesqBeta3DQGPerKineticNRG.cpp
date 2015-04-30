/** 
 * @file BoussinesqBeta3DQGPerKineticNRG.cpp
 * @brief Source of the implementation of the kinetic energy caluclation in the periodic Beta 3DQG model
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBeta3DQGPerKineticNRG.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBeta3DQGPerKineticNRG::BoussinesqBeta3DQGPerKineticNRG(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBeta3DQGPerKineticNRG::~BoussinesqBeta3DQGPerKineticNRG()
   {
   }

   void BoussinesqBeta3DQGPerKineticNRG::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false);
   }

   void BoussinesqBeta3DQGPerKineticNRG::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation:
      ///   \f$ u \cdot u\f$
      ///
      rNLComp.setData((this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::THREE).data().array()*this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::THREE).data().array()).matrix());
      rNLComp.addData((this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::TWO).data().array()*this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::TWO).data().array()).matrix());
      rNLComp.addData((this->scalar(PhysicalNames::VELOCITYZ).dom(0).phys().data().array()*this->scalar(PhysicalNames::VELOCITYZ).dom(0).phys().data().array()).matrix());
   }

   void BoussinesqBeta3DQGPerKineticNRG::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::KINETIC_ENERGY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::KINETIC_ENERGY, FieldRequirement(true, true, false, false));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, false));
   }

}
}
