/** 
 * @file DissB.cpp
 * @brief Source of the implementation of the local ohmic dissipation computation in the F-plane 3DQG model
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/DissB.hpp )

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

namespace QGmhdBhhLowRm {

   DissB::DissB(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   DissB::~DissB()
   {
   }

   void DissB::setCoupling()
   {	
      // 1: want index to start at 1 because of inverse laplacian, T, T?
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false);
   }

   void DissB::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);


      /// 
      /// Computation:
      ///   \f$ FJZ^2 + (\nabla FBZ)^2  \f$
      ///  or
      ///   \f$ (\partial_x b_z)^2 + (\partial_y b_z)^2  + (\partial_y bx - \partial_x by)^2  \f$
      ///
      
      rNLComp.setData((pow(this->scalar(PhysicalNames::FJZ).dom(0).phys().data().array(),2) + pow(this->scalar(PhysicalNames::FBZ).dom(0).grad().comp(FieldComponents::Physical::X).data().array(),2) + pow(this->scalar(PhysicalNames::FBZ).dom(0).grad().comp(FieldComponents::Physical::Y).data().array(),2)).matrix());

   }

   void DissB::setRequirements()
   {
      // Set fluctuating bx field as equation unknown
      this->setName(PhysicalNames::DISSB);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add DissB requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::DISSB, FieldRequirement(true, true, true, false, false, false));

      // Add FJZ requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::FJZ, FieldRequirement(true, true, true, false, false, false));


      // Add FBZ requirements: is scalar?, need spectral?, need physical?, need diff? need curl? need diff2?
      this->mRequirements.addField(PhysicalNames::FBZ, FieldRequirement(true, true, true, true, false, false));

   }
}
}
}
}
}
