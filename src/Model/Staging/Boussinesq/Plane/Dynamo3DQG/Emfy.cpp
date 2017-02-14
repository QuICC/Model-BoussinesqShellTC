/** 
 * @file Emfy.cpp
 * @brief Source of the implementation of EMF y in the Dynamo 3DQG model
 * @author Meredith / Philippe Marti \<philippe.marti@colorado.edu\>
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/Emfy.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace Dynamo3DQG {

   Emfy::Emfy(SharedEquationParameters spEqParams)
      : IScalarTimeAveragedEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Emfy::~Emfy()
   {
   }

   void Emfy::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false);
   }

   void Emfy::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get paramters
//      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);

      /// 
      /// Computation:
      ///   \f$ E_Y \f$
      ///
      //
      rNLComp.setData((this->scalar(PhysicalNames::VELOCITYZ).dom(0).phys().data().array()*this->scalar(PhysicalNames::FBX).dom(0).phys().data().array() + this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::Y).data().array()*this->scalar(PhysicalNames::FBZ).dom(0).phys().data().array()).matrix());

   }

   Datatypes::SpectralScalarType::PointType Emfy::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
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

   void Emfy::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::EMFY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::BEFORE);

      // Add mean temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::EMFY, FieldRequirement(true, true, true, false));

      // Add temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, false, true, false));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, true, true));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::FBZ, FieldRequirement(true, false, true, false));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::FBX, FieldRequirement(true, false, true, false));
   }

}
}
}
}
}
