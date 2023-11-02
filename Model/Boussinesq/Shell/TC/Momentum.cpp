/**
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in
 * the Boussinesq thermal convection in a spherical shell model
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Shell/TC/Momentum.hpp"
#include "Model/Boussinesq/Shell/TC/MomentumKernel.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Transform/Path/I2CurlNl.hpp"
#include "QuICC/Transform/Path/NegI2CurlCurlNl.hpp"
#include "QuICC/Transform/Path/NegI4CurlCurlNl.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace TC {

Momentum::Momentum(SharedEquationParameters spEqParams,
   SpatialScheme::SharedCISpatialScheme spScheme,
   std::shared_ptr<Model::IModelBackend> spBackend) :
    IVectorEquation(spEqParams, spScheme, spBackend)
{
   // Set the variable requirements
   this->setRequirements();
}

Momentum::~Momentum() {}

void Momentum::setCoupling()
{
   int start;
   if (this->ss().has(SpatialScheme::Feature::SpectralOrdering132))
   {
      start = 1;
   }
   else if (this->ss().has(SpatialScheme::Feature::SpectralOrdering123))
   {
      start = 0;
   }
   else
   {
      throw std::logic_error(
         "Unknown spatial scheme was used to setup equations!");
   }

   auto features = defaultCouplingFeature();
   features.at(CouplingFeature::Nonlinear) = true;

   this->defineCoupling(FieldComponents::Spectral::TOR,
      CouplingInformation::PROGNOSTIC, start, features);

   this->defineCoupling(FieldComponents::Spectral::POL,
      CouplingInformation::PROGNOSTIC, start, features);
}

void Momentum::setNLComponents()
{
   this->addNLComponent(FieldComponents::Spectral::TOR,
      Transform::Path::I2CurlNl::id());

   if (this->couplingInfo(FieldComponents::Spectral::POL).isSplitEquation())
   {
      this->addNLComponent(FieldComponents::Spectral::POL,
         Transform::Path::NegI2CurlCurlNl::id());
   }
   else
   {
      this->addNLComponent(FieldComponents::Spectral::POL,
         Transform::Path::NegI4CurlCurlNl::id());
   }
}

void Momentum::initNLKernel(const bool force)
{
   // Initialize if empty or forced
   if (force || !this->mspNLKernel)
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
   this->setSolveTiming(SolveTiming::Prognostic::id());

   // Forward transform generates nonlinear RHS
   this->setForwardPathsType(FWD_IS_NONLINEAR);

   // Get reference to spatial scheme
   const auto& ss = this->ss();

   // Add velocity to requirements: is scalar?
   auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(),
      FieldRequirement(false, ss.spectral(), ss.physical()));
   velReq.enableSpectral();
   velReq.enablePhysical();
   velReq.enableCurl();
}

} // namespace TC
} // namespace Shell
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC
