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
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RandomScalarState::RandomScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mMin(-10), mMax(10), mRatio1D(1e3), mRatio2D(1e3), mRatio3D(1e3)
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

   void RandomScalarState::setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D)
   {
      if(max <= min || ratio1D < 1 || ratio2D < 1 || ratio3D < 1)
      {
         throw Exception("Incompatible spectrum properties requested!");
      }

      // Set min and max range for values
      this->mMin = min;
      this->mMax = max;

      // Set spectrum ratios
      this->mRatio1D = ratio1D;
      this->mRatio2D = ratio2D;
      this->mRatio3D = ratio3D;
   }

   void RandomScalarState::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, false, false, true, false);
      this->setExplicitTiming(FieldComponents::Spectral::SCALAR, ExplicitTiming::LINEAR);
   }

   Datatypes::SpectralScalarType::PointType RandomScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i1D, const int i3D, const int i2D) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      // Get first dimension
      int n1D = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get second dimension
      int n2D = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get third dimension
      int n3D = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      if(i1D < n1D-4 && i3D < n3D - 4 && i2D < n2D - 4)
      {
         // Compute scaling factors
         MHDFloat a1D = exp(-static_cast<MHDFloat>(i1D)*log(this->mRatio1D)/static_cast<MHDFloat>(n1D));
         MHDFloat a2D = exp(-static_cast<MHDFloat>(i2D)*log(this->mRatio2D)/static_cast<MHDFloat>(n2D));
         MHDFloat a3D = exp(-static_cast<MHDFloat>(i3D)*log(this->mRatio3D)/static_cast<MHDFloat>(n3D));

         Datatypes::SpectralScalarType::PointType val;

         this->makeRandom(val, i1D, i3D, i2D);

         return val*a1D*a2D*a3D;
      } else
      {
         return Datatypes::SpectralScalarType::PointType(0);
      }
   }

   void RandomScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

   void RandomScalarState::makeRandom(MHDFloat& val, const int i1D, const int i3D, const int i2D) const
   {
      val = ((this->mMin-this->mMax)*static_cast<MHDFloat>(rand())/RAND_MAX)+this->mMax;
   }

   void RandomScalarState::makeRandom(MHDComplex& val, const int i1D, const int i3D, const int i2D) const
   {
      val.real() = ((this->mMin-this->mMax)*static_cast<MHDFloat>(rand())/RAND_MAX)+this->mMax;

      if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i2D) != 0)
      {
         val.imag() = ((this->mMin-this->mMax)*static_cast<MHDFloat>(rand())/RAND_MAX)+this->mMax;
      } else
      {
         val.imag() = 0.0;
      }
   }

}
}
