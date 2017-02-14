/** 
 * @file TimestepperTFFTimeAvgScalar.cpp
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
#include "Equations/Tests/TimestepperTFFTimeAvgScalar.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"

namespace QuICC {

namespace Equations {

   TimestepperTFFTimeAvgScalar::TimestepperTFFTimeAvgScalar(SharedEquationParameters spEqParams)
      : IScalarTimeAveragedEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   TimestepperTFFTimeAvgScalar::~TimestepperTFFTimeAvgScalar()
   {
   }

   void TimestepperTFFTimeAvgScalar::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, false, false);
   }

   void TimestepperTFFTimeAvgScalar::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      MHDFloat delta = this->eqParams().nd(NonDimensional::DELTA);
      MHDFloat epsilon = this->eqParams().nd(NonDimensional::EPSILON);

      int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);
      Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nJ);
      Array gI = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nI);

      MHDFloat z;
      MHDFloat j_;
      MHDFloat i_;
      nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
      for(int iK = 0; iK < nK; ++iK)
      {
         z = gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
         nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            j_ = gJ(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
            for(int iI = 0; iI < nI; ++iI)
            {
               i_ = gI(iI);

               MHDFloat val = (1./8.)*(std::cos(Math::PI*z/2.0)*(8.0*std::cos(2.0*Math::PI*this->time()) + Math::PI*epsilon*(-1.0 + 2.0*Math::PI + std::sin(2.0*Math::PI*this->time()))) - 2.0*delta*(-1.0 + 2.0*Math::PI + std::sin(2.0*Math::PI*this->time()))*std::sin(Math::PI*z/2.0)) ;
               val -= delta*this->unknown().dom(0).grad().comp(FieldComponents::Physical::Z).point(iI,iJ,iK);
               rNLComp.setPoint(val, iI, iJ, iK);
            }
         }
      }
   }

   void TimestepperTFFTimeAvgScalar::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::PRESSURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::BEFORE);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?, need curl?, need grad2?
      this->mRequirements.addField(PhysicalNames::PRESSURE, FieldRequirement(true, true, true, true, false, false));
   }

}
}
