/** 
 * @file VelocityX.cpp
 * @brief Source of the implementation of the vertical vorticity computation in the F-plane 3DQG model
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/VelocityX.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace QGmhdBhhLowRm {

   VelocityX::VelocityX(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   VelocityX::~VelocityX()
   {
   }

   void VelocityX::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, true, true);
   }

   void VelocityX::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation:
      ///   \f$ VelocityX = -D(streamfunction)Dy \f$
      ///

      rNLComp.setData(-(this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::Y).data().array()).matrix());
   }

   Datatypes::SpectralScalarType::PointType VelocityX::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
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

   void VelocityX::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VELOCITYX);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Set non orthogonal vertical vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYX, FieldRequirement(true, true, true, false));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff? need curl? need diff2?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, true, false, false, true));
   }

}
}
}
}
}
