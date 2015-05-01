/** 
 * @file TestNonlinearScalar.cpp
 * @brief Source of the implementation of a scalar nonlinear test equation
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
#include "Equations/Tests/TestNonlinearScalar.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"

#include "TypeSelectors/TransformSelector.hpp"
namespace GeoMHDiSCC {

namespace Equations {

   TestNonlinearScalar::TestNonlinearScalar(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   TestNonlinearScalar::~TestNonlinearScalar()
   {
   }

   void TestNonlinearScalar::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestNonlinearScalar::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);
      Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nJ);
      Array gI = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nI);

      MHDFloat k_;
      MHDFloat j_;
      MHDFloat i_;
      nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iK = 0; iK < nK; ++iK)
      {
         k_ = gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iK));
         nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            j_ = gJ(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iJ, iK));
            for(int iI = 0; iI < nI; ++iI)
            {
               i_ = gI(iI);

               MHDFloat valK = 1.0;
               MHDFloat valJ = 100*(static_cast<int>(this->name()) + 1)*std::cos(j_);
               MHDFloat valI = 1.0;

               rNLComp.setPoint(valK*valJ*valI, iI, iJ, iK);
            }
         }
      }
   }

   Datatypes::SpectralScalarType::PointType TestNonlinearScalar::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

//      if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iY) == 0)
//      {
//         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iZ,iY) == 0)
//         {
//            if(iX == 0)
//            {
//               return Datatypes::SpectralScalarType::PointType(-1.0);
//            }
//         }
//      } 

      return Datatypes::SpectralScalarType::PointType(0.0);
   }

   void TestNonlinearScalar::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void TestNonlinearScalar::setRequirements()
   {
      // Add requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

}
}
