/** 
 * @file Induction.cpp
 * @brief Source of the implementation of the vector induction equation in the Boussinesq thermal convection dynamo in a sphere model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Sphere/Dynamo/Induction.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   Induction::Induction(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Induction::~Induction()
   {
   }

   void Induction::setCoupling()
   {
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         int start = 1;
      #elif defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         int start = 0;
      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);
   }

   void Induction::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::POL,0);

      this->addNLComponent(FieldComponents::Spectral::TOR,1);
   }

   void Induction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      ///
      /// Compute \f$\left(\vec u\wedge\vec B\right)\f$
      ///
      switch(compId)
      {
         case(FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), 1.0);
            break;
         default:
            assert(false);
            break;
      }
   }

   void Induction::setRequirements()
   {
      // Set velocity as equation unknown
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
}
}
}