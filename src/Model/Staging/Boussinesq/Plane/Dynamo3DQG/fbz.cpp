/** 
 * @file fbz.cpp
 * @brief Source of the implementation of the mean heat computation in the F-plane 3DQG model
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/Dynamo3DQG/fbz.hpp )

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

   fbz::fbz(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   fbz::~fbz()
   {
   }

   void fbz::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, true, true);
   }

   void fbz::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get paramters
      MHDFloat MPr = this->eqParams().nd(NonDimensional::MAGPRANDTL);

      /// 
      /// Computation:
      ///
      rNLComp.setData((-MPr*(this->scalar(PhysicalNames::BX).dom(0).phys().data().array()*this->scalar(PhysicalNames::VELOCITYZ).dom(0).grad().comp(FieldComponents::Physical::X).data().array()+this->scalar(PhysicalNames::BY).dom(0).phys().data().array()*this->scalar(PhysicalNames::VELOCITYZ).dom(0).grad().comp(FieldComponents::Physical::Y).data().array())).matrix());
   }

   Datatypes::SpectralScalarType::PointType fbz::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
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

   void fbz::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::FBZ);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add fbz requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::FBZ, FieldRequirement(true, true, true, false));

      // Add BX requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::BX, FieldRequirement(true, true, true, false));

      // Add BY requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::BY, FieldRequirement(true, true, true, false));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));
   }

}
}
}
}
}
