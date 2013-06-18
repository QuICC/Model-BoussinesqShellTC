/** \file RandomScalarState.cpp
 *  \brief Source of the implementation of the general random scalar state equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Generator/States/RandomScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RandomScalarState::RandomScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mMin(-10), mMax(10), mXRatio(1e3), mYRatio(1e3), mZRatio(1e3)
   {
   }

   RandomScalarState::~RandomScalarState()
   {
   }

   void RandomScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void RandomScalarState::setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat xRatio, const MHDFloat yRatio, const MHDFloat zRatio)
   {
      if(max <= min || xRatio < 1 || yRatio < 1 || zRatio < 1)
      {
         throw Exception("Incompatible spectrum properties requested!");
      }

      // Set min and max range for values
      this->mMin = min;
      this->mMax = max;

      // Set spectrum ratios
      this->mXRatio = xRatio;
      this->mYRatio = yRatio;
      this->mZRatio = zRatio;
   }

   void RandomScalarState::setCoupling()
   {
      // Get X dimension
      int nX = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get Y dimension
      int nY = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: first complex solver, complex solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, false, 0);

      // Set nonlinear flags: NO nonlinear term, NO quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: has source term
      infoIt.first->second.setSource(true);

      // Equation is coupled to itself
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Set sizes of blocks and matrices
      ArrayI blockNs(nY);
      blockNs.setConstant(nX*nZ);
      ArrayI rhsCols(nY);
      rhsCols.setConstant(1);
      infoIt.first->second.setSizes(nY, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   MHDComplex RandomScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      // Get X dimension
      int nX = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get Y dimension
      int nY = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      if(iX < nX-4 && iZ < nZ - 4 && iY < nY - 4)
      {
         // Compute scaling factors
         MHDFloat aX = exp(-static_cast<MHDFloat>(iX)*log(this->mXRatio)/static_cast<MHDFloat>(nX));
         MHDFloat aY = exp(-static_cast<MHDFloat>(iY)*log(this->mYRatio)/static_cast<MHDFloat>(nY));
         MHDFloat aZ = exp(-static_cast<MHDFloat>(iZ)*log(this->mZRatio)/static_cast<MHDFloat>(nZ));

         MHDComplex val;
         val.real() = ((this->mMin-this->mMax)*static_cast<MHDFloat>(rand())/RAND_MAX)+this->mMax;
         val.imag() = ((this->mMin-this->mMax)*static_cast<MHDFloat>(rand())/RAND_MAX)+this->mMax;

         return val*aX*aY*aZ;
      } else
      {
         return MHDComplex(0,0);
      }
   }

   void RandomScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

}
}
