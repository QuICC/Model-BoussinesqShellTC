/** 
 * @file MeanHeat.cpp
 * @brief Source of the implementation of the mean heat equation in the periodic Beta 3DQG model
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Beta3DQG/MeanHeat.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "PhysicalOperators/StreamHeatAdvection.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace Beta3DQG {

   MeanHeat::MeanHeat(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   MeanHeat::~MeanHeat()
   {
   }

   void MeanHeat::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false);
   }

   void MeanHeat::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get paramters
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);

      /// 
      /// Computation:
      ///   \f$ Pr u \theta\f$
      ///
      rNLComp.setData((-Pr*this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::Y).data().array()*this->scalar(PhysicalNames::TEMPERATURE).dom(0).phys().data().array()).matrix());
   }

   void MeanHeat::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::DX_MEANTEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::DX_MEANTEMPERATURE, FieldRequirement(true, true, false, false));

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, false, true, false));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));

      // Gradient does not require Z component
      ArrayB   comps = ArrayB::Constant(3, true);
      comps(0) = false;
      std::map<FieldComponents::Spectral::Id,ArrayB>  gradComps;
      gradComps.insert(std::make_pair(FieldComponents::Spectral::SCALAR, comps));

      // Update streamfunction gradient requirements
      this->updateFieldRequirements(PhysicalNames::STREAMFUNCTION).updateGradient(gradComps);
   }

}
}
}
}
}