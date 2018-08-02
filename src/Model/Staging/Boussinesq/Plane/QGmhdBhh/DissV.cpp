/** 
 * @file DissV.cpp
 * @brief Source of the implementation of the local viscous dissipation computation in the F-plane 3DQG model
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

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhh/DissV.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "Enums/FieldIdsTools.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace QGmhdBhh {

   DissV::DissV(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   DissV::~DissV()
   {
   }

   void DissV::setCoupling()
   {	
      // 1: want index to start at 1 because of inverse laplacian, T, T?
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false);
   }

   void DissV::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);


      /// 
      /// Computation:
      ///   \f$ VORTICITYZ^2 + (\nabla VELOCITYZ)^2  \f$
      ///

        rNLComp.setData((pow(this->scalar(PhysicalNames::VORTICITYZ).dom(0).phys().data().array(),2) + pow(this->scalar(PhysicalNames::VELOCITYZ).dom(0).grad().comp(FieldComponents::Physical::X).data().array(),2) + pow(this->scalar(PhysicalNames::VELOCITYZ).dom(0).grad().comp(FieldComponents::Physical::Y).data().array(),2)).matrix());

   }

   void DissV::setRequirements()
   {
      // Set fluctuating bx field as equation unknown
      this->setName(PhysicalNames::DISSV);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add DissV requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::DISSV, FieldRequirement(true, true, true, false, false, false));

      // Add VORTICITYZ requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, true, false, false, false));


      // Add VELOCITYZ requirements: is scalar?, need spectral?, need physical?, need diff? need curl? need diff2?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true, false, false));

   }
}
}
}
}
}
