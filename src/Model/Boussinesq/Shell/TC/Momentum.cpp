/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq thermal convection in a spherical shell model
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Shell/TC/Momentum.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SpatialScheme/3D/SLFl.hpp"
#include "QuICC/SpatialScheme/3D/SLFm.hpp"
#include "QuICC/Model/Boussinesq/Shell/TC/MomentumKernel.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace TC {

   Momentum::Momentum(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IVectorEquation(spEqParams,spScheme)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Momentum::~Momentum()
   {
   }

   void Momentum::setCoupling()
   {
      int start;
      if(this->ss().id() == SpatialScheme::SLFl::sId)
      {
         start = 1;
      } else if(this->ss().id() == SpatialScheme::SLFm::sId)
      {
         start = 0;
      } else
      {
         throw std::logic_error("Unknown spatial scheme was used to setup equations!");
      }

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);
   }

   void Momentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, 0);

      this->addNLComponent(FieldComponents::Spectral::POL, 0);
   }

   void Momentum::initNLKernel(const bool force)
   {
      // Initialize if empty or forced
      if(force || !this->mspNLKernel)
      {
         // Initialize the physical kernel
         auto spNLKernel = std::make_shared<Physical::Kernel::MomentumKernel>();
         spNLKernel->setVelocity(this->name(), this->spUnknown());
         spNLKernel->init(1.0);
         this->mspNLKernel = spNLKernel;
      }
   }

   void Momentum::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::Velocity::id());

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Forward transform generates nonlinear RHS
      this->setForwardPathsType(FWD_IS_NONLINEAR);

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add velocity to requirements: is scalar?
      auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(), FieldRequirement(false, ss.spectral(), ss.physical()));
      velReq.enableSpectral();
      velReq.enablePhysical();
      velReq.enableCurl();
   }

}
}
}
}
}
