/** 
 * @file BoussinesqRBCSquareSMomentum.cpp
 * @brief Source of the implementation of the streamfunction momentum equation in Rayleigh-Benard convection in a square cavity (2D) (streamfunction formulation)
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
#include "Equations/Box/Boussinesq/BoussinesqRBCSquareSMomentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"

namespace QuICC {

namespace Equations {

   BoussinesqRBCSquareSMomentum::BoussinesqRBCSquareSMomentum(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqRBCSquareSMomentum::~BoussinesqRBCSquareSMomentum()
   {
   }

   void BoussinesqRBCSquareSMomentum::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void BoussinesqRBCSquareSMomentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Get Prandtl number
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)u_x\f$
      ///
      Physical::StreamAdvection<FieldComponents::Physical::X,FieldComponents::Physical::Z>::set(rNLComp, this->unknown().dom(0).grad(), this->scalar(PhysicalNames::VORTICITYY).dom(0).grad(), 1.0/Pr);
   }

   void BoussinesqRBCSquareSMomentum::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::STREAMFUNCTION);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, true));

      // Add vorticityz to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VORTICITYY, FieldRequirement(true, true, false, true));
   }

}
}
