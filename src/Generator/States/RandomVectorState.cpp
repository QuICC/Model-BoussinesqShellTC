/** 
 * @file RandomVectorState.cpp
 * @brief Source of the implementation of the general random scalar state equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <time.h>

// External includes
//

// Class include
//
#include "Generator/States/RandomVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"

namespace QuICC {

namespace Equations {

   RandomVectorState::RandomVectorState(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mStartSeed(1)
   {
      timespec t;
      clock_gettime(CLOCK_REALTIME, &t);
      srand(t.tv_sec + t.tv_nsec);

      this->mStartSeed = t.tv_sec + t.tv_nsec;
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

   void RandomVectorState::setSpectrum(const FieldComponents::Spectral::Id compId, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D)
   {
      if(max <= min || ratio1D < 1 || ratio2D < 1)
      {
         throw Exception("Incompatible spectrum properties requested!");
      }

      // Set min and max range for values
      this->mMin.insert(std::make_pair(compId, min));
      this->mMax.insert(std::make_pair(compId, max));

      // Set spectrum ratios
      this->mRatio1D.insert(std::make_pair(compId, ratio1D));
      this->mRatio2D.insert(std::make_pair(compId, ratio2D));
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
         this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::TRIVIAL, 0, false, true, false, false);
      }
      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::TRIVIAL, 0, false, true, false, false);
      }
      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::THREE, CouplingInformation::TRIVIAL, 0, false, true, false, false);
      }
   }

   Datatypes::SpectralScalarType::PointType RandomVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      MHDFloat minVal = this->mMin.find(compId)->second;
      MHDFloat maxVal = this->mMax.find(compId)->second;
      MHDFloat ratio1D = this->mRatio1D.find(compId)->second;
      MHDFloat ratio2D = this->mRatio2D.find(compId)->second;

      #ifdef QUICC_SPATIALDIMENSION_3D
         // Get ratios for components
         MHDFloat ratio3D = this->mRatio3D.find(compId)->second;

         // Get first dimension
         int n1D = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         // Get second dimension
         int n2D = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
         // Get third dimension
         int n3D = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

         // Get simulation wide indexes
         int j_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
         int k_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

         int z1D = 4;
         int z2D = 4;
         int z3D = 4;
         #if defined QUICC_SPATIALSCHEME_TFF
         z2D = 2;
         if(k_ >= n2D/2)
         {
            k_ = n2D - k_;
         }
         n2D = n2D/2;
         #endif //defined QUICC_SPATIALSCHEME_TFF

         if(i < n1D-z1D && j_ < n3D - z3D && k_ < n2D - z2D)
         {
            // Compute scaling factors
            MHDFloat a1D = exp(-static_cast<MHDFloat>(i)*log(ratio1D)/static_cast<MHDFloat>(n1D));
            MHDFloat a2D = exp(-static_cast<MHDFloat>(k_)*log(ratio2D)/static_cast<MHDFloat>(n2D));
            MHDFloat a3D = exp(-static_cast<MHDFloat>(j_)*log(ratio3D)/static_cast<MHDFloat>(n3D));

            Datatypes::SpectralScalarType::PointType val;

            this->makeRandom(val, i, j, k, minVal, maxVal);

            return val*a1D*a2D*a3D;
         } else
         {
            return Datatypes::SpectralScalarType::PointType(0);
         }
      #else
         // Get first dimension
         int n1D = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         // Get second dimension
         int n2D = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

         // Get simulation wide indexes
         int j_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j);

         int z1D = 4;
         int z2D = 4;

         if(i < n1D-z1D && j_ < n2D - z2D)
         {
            // Compute scaling factors
            MHDFloat a1D = exp(-static_cast<MHDFloat>(i)*log(ratio1D)/static_cast<MHDFloat>(n1D));
            MHDFloat a2D = exp(-static_cast<MHDFloat>(j_)*log(ratio2D)/static_cast<MHDFloat>(n2D));

            Datatypes::SpectralScalarType::PointType val;

            this->makeRandom(val, i, j, k, minVal, maxVal);

            return val*a1D*a2D;
         } else
         {
            return Datatypes::SpectralScalarType::PointType(0);
         }
      #endif //QUICC_SPATIALDIMENSION_3D
   }

   void RandomVectorState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, true, false));
   }

   void RandomVectorState::makeRandom(MHDFloat& val, const int i, const int j, const int k, const MHDFloat minVal, const MHDFloat maxVal, const unsigned int seed) const
   {
      if(seed == 1)
      {
         val = ((minVal-maxVal)*static_cast<MHDFloat>(rand())/RAND_MAX)+maxVal;
      } else
      {
         srand(this->mStartSeed + seed);
         val = ((minVal-maxVal)*static_cast<MHDFloat>(rand())/RAND_MAX)+maxVal;
         timespec t;
         clock_gettime(CLOCK_REALTIME, &t);
         srand(t.tv_sec + t.tv_nsec);
      }
   }

   void RandomVectorState::makeRandom(MHDComplex& val, const int i, const int j, const int k, const MHDFloat minVal, const MHDFloat maxVal) const
   {
      MHDFloat tmp;
      this->makeRandom(tmp, i, j, k, minVal, maxVal);

      val.real() = tmp;

      #if defined QUICC_SPATIALSCHEME_TFT || defined QUICC_SPATIALSCHEME_FFF || defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM || defined QUICC_SPATIALSCHEME_AFT || defined QUICC_SPATIALSCHEME_CFT || defined QUICC_SPATIALSCHEME_WFT
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) != 0)
         {
            this->makeRandom(tmp, i, j, k, minVal, maxVal);
            val.imag() = tmp;
         } else
         {
            val.imag() = 0.0;
         }
      #elif defined QUICC_SPATIALSCHEME_TFF
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 0 && this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) != 0)
         {
            unsigned int seed = 2;
            seed += this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATF1D>(i,k);

            int n2D = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
            int k2D = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if(k2D < n2D/2)
            {
               seed += k2D;
            } else
            {
               seed += n2D - k2D;
            }

            this->makeRandom(tmp, i, j, k, minVal, maxVal, seed);
            val.real() = tmp;
            this->makeRandom(tmp, i, j, k, minVal, maxVal, seed+3);
            if(k2D < n2D/2)
            {
               val.imag() = tmp;
            } else
            {
               val.imag() = -tmp;
            }

         } else if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) != 0)
         {
            this->makeRandom(tmp, i, j, k, minVal, maxVal);
            val.imag() = tmp;
         } else
         {
            val.imag() = 0.0;
         }
      #elif defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k) != 0)
         {
            this->makeRandom(tmp, i, j, k, minVal, maxVal);
            val.imag() = tmp;
         } else
         {
            val.imag() = 0.0;
         }
      #endif //defined QUICC_SPATIALSCHEME_TFT || defined QUICC_SPATIALSCHEME_FFF || defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM || defined QUICC_SPATIALSCHEME_AFT || defined QUICC_SPATIALSCHEME_CFT || defined QUICC_SPATIALSCHEME_WFT
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
