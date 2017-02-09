/** 
 * @file BoussinesqDynamoPlaneInduction.cpp
 * @brief Source of the implementation of the vector induction equation for rotating thermal dynamo in a plane layer (toroidal/poloidal formulation)
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
#include "Equations/Box/Boussinesq/BoussinesqDynamoPlaneInduction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Equations {

   BoussinesqDynamoPlaneInduction::BoussinesqDynamoPlaneInduction(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqDynamoPlaneInduction::~BoussinesqDynamoPlaneInduction()
   {
   }

   void BoussinesqDynamoPlaneInduction::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, 0, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqDynamoPlaneInduction::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::POL, 1);

      this->addNLComponent(FieldComponents::Spectral::TOR, 1);
   }

   void BoussinesqDynamoPlaneInduction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      ///
      /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
      ///
      switch(compId)
      {
         case(FieldComponents::Physical::X):
            Physical::Cross<FieldComponents::Physical::Y,FieldComponents::Physical::Z>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::Y):
            Physical::Cross<FieldComponents::Physical::Z,FieldComponents::Physical::X>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::Z):
            Physical::Cross<FieldComponents::Physical::X,FieldComponents::Physical::Y>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), 1.0);
            break;
         default:
            assert(false);
            break;
      }
   }

   void BoussinesqDynamoPlaneInduction::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::MAGNETIC);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::MAGNETIC, FieldRequirement(false, true, true, false, false));

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, false));
   }

}
}
