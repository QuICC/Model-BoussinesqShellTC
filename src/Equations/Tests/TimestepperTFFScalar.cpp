/** 
 * @file TimestepperTFFScalar.cpp
 * @brief Source of the implementation of the equation for the timestepper in TFF scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "TypeSelectors/TransformSelector.hpp"

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Tests/TimestepperTFFScalar.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TimestepperTFFScalar::TimestepperTFFScalar(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   TimestepperTFFScalar::~TimestepperTFFScalar()
   {
   }

   void TimestepperTFFScalar::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void TimestepperTFFScalar::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      MHDFloat eps = this->eqParams().nd(NonDimensional::EPSILON);

      int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);
      Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nJ);
      Array gI = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nI);

      MHDFloat k_;
      MHDFloat j_;
      MHDFloat i_;
      MHDFloat z;
      nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
      for(int iK = 0; iK < nK; ++iK)
      {
         k_ = gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
         z = (k_ + 1.0)/2.0;
         nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            j_ = gJ(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
            for(int iI = 0; iI < nI; ++iI)
            {
               i_ = gI(iI);

               MHDFloat val = 2.0*(z - 1.0)*z*(147.0 + z*(-714.0 + z*(589.0 + 10.0*z*(16.0 + z))))*eps/std::pow(7 + z,3);

               rNLComp.setPoint(val, iI, iJ, iK);
            }
         }
      }

   }

   void TimestepperTFFScalar::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?, need curl?, need grad2?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, false));
   }

}
}
