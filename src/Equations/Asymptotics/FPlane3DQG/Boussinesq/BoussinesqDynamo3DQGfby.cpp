/** 
 * @file BoussinesqDynamo3DQGfby.cpp
 * @brief Source of the implementation of the mean heat computation in the F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGfby.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqDynamo3DQGfby::BoussinesqDynamo3DQGfby(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqDynamo3DQGfby::~BoussinesqDynamo3DQGfby()
   {
   }

   void BoussinesqDynamo3DQGfby::setCoupling()
   {	
      // 1: want index to start at 1 because of inverse laplacian, T, T?
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, true, true);
   }

   void BoussinesqDynamo3DQGfby::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get paramters
      MHDFloat MPr = this->eqParams().nd(NonDimensional::MAGPRANDTL);

      /// 
      /// Computation:
      ///   \f$ MPr(Bx dxy\Psi + By dyy\Psi) \f$
      ///
      rNLComp.setData((-MPr*(this->scalar(PhysicalNames::BX).dom(0).phys().data().array()*this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad().comp(FieldComponents::Physical::X,X).data().array() + this->scalar(PhysicalNames::BY).dom(0).phys().data().array*this->scalar(PhysicalNames:STREAMFUNCTION).dom().grad().comp(FieldComponents::Physical::X,Y).data().array())).matrix());
   }

   Datatypes::SpectralScalarType::PointType BoussinesqDynamo3DQGfby::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
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

   void BoussinesqDynamo3DQGfby::setRequirements()
   {
      // Set fluctuating bx field as equation unknown
      this->setName(PhysicalNames::FBY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add fBY requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::FBY, FieldRequirement(true, true, true, false));

      // Add BX requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::BX, FieldRequirement(true, true, true, false));

      // Add BY requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::BY, FieldRequirement(true, true, true, false));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, true, false));
   }

}
}
