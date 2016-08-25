/** 
 * @file BoussinesqDynamo3DQGEmfx.cpp
 * @brief Source of the implementation of the mean Bx in the Dynamo 3DQG model
 * @author Meredith / Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGEmfx.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqDynamo3DQGEmfx::BoussinesqDynamo3DQGBx(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqDynamo3DQGEmfx::~BoussinesqDynamo3DQGBx()
   {
   }

   void BoussinesqDynamo3DQGEmfx::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqDynamo3DQGEmfx::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get paramters
//      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);

      /// 
      /// Computation:
      ///   \f$  EMFx \f$
      ///

      rNLComp.setData((-this->scalar(PhysicalNames::VELOCITYZ).dom(0).phys().data().array()*this->scalar(PhysicalNames::FBY).dom(0).phys().data().array() + this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::X).data().array()*this->scalar(PhysicalNames::FBZ).dom(0).phys().data().array()).matrix());

   }

   Datatypes::SpectralScalarType::PointType BoussinesqDynamo3DQGEmfx::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
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

   void BoussinesqDynamo3DQGEmfx::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::EMFX);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add mean temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::EMFX, FieldRequirement(true, true, true, false));

      // Add temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, false, true, false));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, true, true));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::FBZ, FieldRequirement(true, false, true, false));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::FBY, FieldRequirement(true, false, true, false));
   }

}
}
