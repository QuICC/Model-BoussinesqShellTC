/** 
 * @file RandomVectorState.cpp
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
#include "Generator/States/RandomVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RandomVectorState::RandomVectorState(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
   }

   RandomVectorState::~RandomVectorState()
   {
   }

   void RandomVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void RandomVectorState::setSpectrum(const FieldComponents::Spectral::Id compId, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D)
   {
      if(max <= min || ratio1D < 1 || ratio2D < 1 || ratio3D < 1)
      {
         throw Exception("Incompatible spectrum properties requested!");
      }

      // Set min and max range for values
      this->mMin.insert(std::make_pair(compId, min));
      this->mMax.insert(std::make_pair(compId, max));

      // Set spectrum ratios
      this->mRatio1D.insert(std::make_pair(compId, ratio1D));
      this->mRatio2D.insert(std::make_pair(compId, ratio2D));
      this->mRatio3D.insert(std::make_pair(compId, ratio3D));
   }

   void RandomVectorState::setCoupling()
   {
      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::TRIVIAL, 0, false, true, false);
      }
      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::TRIVIAL, 0, false, true, false);
      }
      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::THREE, CouplingInformation::TRIVIAL, 0, false, true, false);
      }
   }

   Datatypes::SpectralScalarType::PointType RandomVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i1D, const int i3D, const int i2D) const
   {
      // Get ratios for components
      MHDFloat minVal = this->mMin.find(compId)->second;
      MHDFloat maxVal = this->mMax.find(compId)->second;
      MHDFloat ratio1D = this->mRatio1D.find(compId)->second;
      MHDFloat ratio2D = this->mRatio2D.find(compId)->second;
      MHDFloat ratio3D = this->mRatio3D.find(compId)->second;

      // Get first dimension
      int n1D = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get second dimension
      int n2D = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get third dimension
      int n3D = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      if(i1D < n1D-4 && i3D < n3D - 4 && i2D < n2D - 4)
      {
         // Compute scaling factors
         MHDFloat a1D = exp(-static_cast<MHDFloat>(i1D)*log(ratio1D)/static_cast<MHDFloat>(n1D));
         MHDFloat a2D = exp(-static_cast<MHDFloat>(i2D)*log(ratio2D)/static_cast<MHDFloat>(n2D));
         MHDFloat a3D = exp(-static_cast<MHDFloat>(i3D)*log(ratio3D)/static_cast<MHDFloat>(n3D));

         Datatypes::SpectralScalarType::PointType val;

         this->makeRandom(val, i1D, i3D, i2D, minVal, maxVal);

         return val*a1D*a2D*a3D;
      } else
      {
         return Datatypes::SpectralScalarType::PointType(0);
      }
   }

   void RandomVectorState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, true, false));
   }

   void RandomVectorState::makeRandom(MHDFloat& val, const int i1D, const int i3D, const int i2D, const MHDFloat minVal, const MHDFloat maxVal) const
   {
      val = ((minVal-maxVal)*static_cast<MHDFloat>(rand())/RAND_MAX)+maxVal;
   }

   void RandomVectorState::makeRandom(MHDComplex& val, const int i1D, const int i3D, const int i2D, const MHDFloat minVal, const MHDFloat maxVal) const
   {
      MHDFloat tmp;
      this->makeRandom(tmp, i1D, i3D, i2D, minVal, maxVal);

      val.real() = tmp;

      #if defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMDHDISCC_SPATIALSCHEME_TFF || defined GEOMDHDISCC_SPATIALSCHEME_FFF || defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM || defined GEOMHDISCC_SPATIALSCHEME_AFT || defined GEOMHDISCC_SPATIALSCHEME_CFT
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i2D) != 0)
         {
            this->makeRandom(tmp, i1D, i3D, i2D, minVal, maxVal);
            val.imag() = tmp;
         } else
         {
            val.imag() = 0.0;
         }
      #elif defined GEOMHDISCC_SPATIALSCHEME_SLFL || defined GEOMHDISCC_SPATIALSCHEME_BLFL
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(i3D, i2D) != 0)
         {
            this->makeRandom(tmp, i1D, i3D, i2D, minVal, maxVal);
            val.imag() = tmp;
         } else
         {
            val.imag() = 0.0;
         }
      #endif //defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMDHDISCC_SPATIALSCHEME_TFF || defined GEOMDHDISCC_SPATIALSCHEME_FFF || defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM || defined GEOMHDISCC_SPATIALSCHEME_AFT || defined GEOMHDISCC_SPATIALSCHEME_CFT
   }

   void RandomVectorState::setNLComponents()
   {
      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(FieldComponents::Spectral::ONE, 0);
      }

      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(FieldComponents::Spectral::TWO, 0);
      }

      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->addNLComponent(FieldComponents::Spectral::THREE, 0);
      }
   }

}
}
