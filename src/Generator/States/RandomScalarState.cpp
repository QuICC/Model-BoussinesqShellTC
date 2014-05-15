/** 
 * @file RandomScalarState.cpp
 * @brief Source of the implementation of the general random scalar state equation
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
#include "Generator/States/RandomScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RandomScalarState::RandomScalarState(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName, spEqParams), mMin(-10), mMax(10), mXRatio(1e3), mYRatio(1e3), mZRatio(1e3)
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
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, false, false, true);
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
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iY) != 0)
         {
            val.imag() = ((this->mMin-this->mMax)*static_cast<MHDFloat>(rand())/RAND_MAX)+this->mMax;
         } else
         {
            val.imag() = 0.0;
         }

         return val*aX*aY*aZ;
      } else
      {
         return MHDComplex(0,0);
      }
   }

   void RandomScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

}
}
