/** \file Beta3DQGTransport.cpp
 *  \brief Source of the implementation of the transport equation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   Beta3DQGTransport::Beta3DQGTransport(SharedIEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Beta3DQGTransport::~Beta3DQGTransport()
   {
   }

   void Beta3DQGTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T}\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);
   }

   void Beta3DQGTransport::computeLinear(Datatypes::SpectralScalarType& rRHS) const
   {  
      ///
      /// Compute the background state advection term
      ///

      // Compute the complex box scale factor
      MHDComplex c = -this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*MathConstants::cI;

      // Get size of dealiased output (at this stage data has still dealiasing rows)
      int dealiasedRows = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Loop over m
      MHDFloat m_;
      int nSlice = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      for(int m = 0; m < nSlice; m++)
      {
         m_ = 0.5*static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m));

         rRHS.addSlice((m_*c)*this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).perturbation().slice(m), m, dealiasedRows);
      }
   }

   void Beta3DQGTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, true));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

   void Beta3DQGTransport::setCoupling()
   {
      // Get X dimension
      int nX = this->unknown().dom(0).spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get Y dimension
      int nY = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Set coupling information
      this->mCouplingInfos.insert(make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      Beta3DQGSystem::setCouplingInfo(this->mCouplingInfos.find(FieldComponents::Spectral::SCALAR)->second, this->name(), nX, nZ, nY);
   }
}
}
