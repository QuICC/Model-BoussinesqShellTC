/** \file Beta3DQGStreamfunction.cpp
 *  \brief Source of the implementation of the streamfunction equation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGStreamfunction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGSystem.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   Beta3DQGStreamfunction::Beta3DQGStreamfunction(SharedIEquationParameters spEqParams)
      : IBeta3DQGScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Beta3DQGStreamfunction::~Beta3DQGStreamfunction()
   {
   }

   void Beta3DQGStreamfunction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\nabla^2_{\perp}\psi\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->unknown().dom(0).grad(), this->scalar(PhysicalNames::VORTICITYZ).dom(0).grad(), 1.0);
   }

   void Beta3DQGStreamfunction::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::STREAMFUNCTION);

      // Set streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, true, true));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, false, true));

      // Add temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, false, false, true));
   }

   void Beta3DQGStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarEquation::timestepOutput(id, storage, matIdx, start);

      // Get the box scale
      MHDFloat boxScale = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);

      // Get right wave number
      MHDFloat m_ = boxScale*0.5*static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx));

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL));

      ///
      /// Compute the vertical vorticity: \f$\zeta = \nabla^2\psi\f$
      ///
      this->rScalar(PhysicalNames::VORTICITYZ).rDom(0).rPerturbation().setSlice(Spectral::PeriodicOperator::laplacian2D(spec1D, m_, 0)*this->unknown().dom(0).perturbation().slice(matIdx), matIdx);
   }

   void Beta3DQGStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarEquation::timestepOutput(id, storage, matIdx, start);

      // Get the box scale
      MHDFloat boxScale = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);

      // Get right wave number
      MHDFloat m_ = boxScale*0.5*static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx));

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL));

      ///
      /// Compute the vertical vorticity: \f$\zeta = \nabla^2\psi\f$
      ///
      this->rScalar(PhysicalNames::VORTICITYZ).rDom(0).rPerturbation().setSlice(Spectral::PeriodicOperator::laplacian2D(spec1D, m_, 0)*this->unknown().dom(0).perturbation().slice(matIdx), matIdx);
   }

}
}
