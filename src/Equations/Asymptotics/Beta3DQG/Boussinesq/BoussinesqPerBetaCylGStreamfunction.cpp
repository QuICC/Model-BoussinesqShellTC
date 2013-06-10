/** \file BoussinesqPerBetaCylGStreamfunction.cpp
 *  \brief Source of the implementation of the streamfunction equation in the 3DQG beta model with periodic radius
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGStreamfunction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGSystem.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqPerBetaCylGStreamfunction::BoussinesqPerBetaCylGStreamfunction(SharedEquationParameters spEqParams)
      : IBoussinesqPerBetaCylGScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqPerBetaCylGStreamfunction::~BoussinesqPerBetaCylGStreamfunction()
   {
   }

   void BoussinesqPerBetaCylGStreamfunction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\nabla^2_{\perp}\psi\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->unknown().dom(0).grad(), this->scalar(PhysicalNames::VORTICITYZ).dom(0).grad(), 0.0);
   }

   void BoussinesqPerBetaCylGStreamfunction::setRequirements()
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

   void BoussinesqPerBetaCylGStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarPEquation::timestepOutput(id, storage, matIdx, start);

      ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get radial wave number rescaled to box size
      MHDFloat kX_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*0.5*static_cast<MHDFloat>(mode(1));
      // Get azimuthal wave number rescaled to box size
      MHDFloat kY_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*0.5*static_cast<MHDFloat>(mode(0));

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL));

      ///
      /// Compute the vertical vorticity: \f$\zeta = \nabla^2\psi\f$
      ///
      this->rScalar(PhysicalNames::VORTICITYZ).rDom(0).rPerturbation().setProfile(Spectral::PeriodicOperator::laplacian2D(kX_, kY_)*spec1D.id(0)*this->unknown().dom(0).perturbation().profile(mode(1),mode(0)), mode(1), mode(0));
   }

   void BoussinesqPerBetaCylGStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarPEquation::timestepOutput(id, storage, matIdx, start);

      ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get radial wave number rescaled to box size
      MHDFloat kX_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*0.5*static_cast<MHDFloat>(mode(1));
      // Get azimuthal wave number rescaled to box size
      MHDFloat kY_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*0.5*static_cast<MHDFloat>(mode(0));

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL));

      ///
      /// Compute the vertical vorticity: \f$\zeta = \nabla^2\psi\f$
      ///
      this->rScalar(PhysicalNames::VORTICITYZ).rDom(0).rPerturbation().setProfile(Spectral::PeriodicOperator::laplacian2D(kX_, kY_)*spec1D.id(0)*this->unknown().dom(0).perturbation().profile(mode(1),mode(0)), mode(1), mode(0));
   }

}
}
