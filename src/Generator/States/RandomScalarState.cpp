/** 
 * @file RandomScalarState.cpp
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
#include "Generator/States/RandomScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RandomScalarState::RandomScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mMin(-10), mMax(10), mRatio1D(1e3), mRatio2D(1e3), mRatio3D(1e3), mStartSeed(1), mSpecial(NOTHING)
   {
      timespec t;
      clock_gettime(CLOCK_REALTIME, &t);
      srand(t.tv_sec + t.tv_nsec);

      this->mStartSeed = t.tv_sec + t.tv_nsec;
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

   void RandomScalarState::setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const SpecialId id)
   {
      if(max <= min || ratio1D < 1 || ratio2D < 1)
      {
         throw Exception("Incompatible spectrum properties requested!");
      }

      // Set min and max range for values
      this->mMin = min;
      this->mMax = max;

      // Set spectrum ratios
      this->mRatio1D = ratio1D;
      this->mRatio2D = ratio2D;

      this->mSpecial = id;
   }

   void RandomScalarState::setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D, const SpecialId id)
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

      this->mSpecial = id;
   }

   void RandomScalarState::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, false, true, false);
   }

   Datatypes::SpectralScalarType::PointType RandomScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      #ifdef QUICC_SPATIALDIMENSION_3D
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
         if(this->mSpecial == ONLYMEAN && !(j_ == 0 && k_ == 0))
         {
            return Datatypes::SpectralScalarType::PointType(0);
         } else if(this->mSpecial == NOMEAN && j_ == 0 && k_ == 0)
         {
            return Datatypes::SpectralScalarType::PointType(0);
         }
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
            MHDFloat a1D = exp(-static_cast<MHDFloat>(i)*log(this->mRatio1D)/static_cast<MHDFloat>(n1D));
            MHDFloat a2D = exp(-static_cast<MHDFloat>(k_)*log(this->mRatio2D)/static_cast<MHDFloat>(n2D));
            MHDFloat a3D = exp(-static_cast<MHDFloat>(j_)*log(this->mRatio3D)/static_cast<MHDFloat>(n3D));

            Datatypes::SpectralScalarType::PointType val;

            this->makeRandom(val, i, j, k);

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
            MHDFloat a1D = exp(-static_cast<MHDFloat>(i)*log(this->mRatio1D)/static_cast<MHDFloat>(n1D));
            MHDFloat a2D = exp(-static_cast<MHDFloat>(j_)*log(this->mRatio2D)/static_cast<MHDFloat>(n2D));

            Datatypes::SpectralScalarType::PointType val;

            this->makeRandom(val, i, j, k);

            return val*a1D*a2D;
         } else
         {
            return Datatypes::SpectralScalarType::PointType(0);
         }
      #endif //QUICC_SPATIALDIMENSION_3D
   }

   void RandomScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::BEFORE);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

   void RandomScalarState::makeRandom(MHDFloat& val, const int i, const int j, const int k, const unsigned int seed) const
   {
      if(seed == 1)
      {
         val = ((this->mMin-this->mMax)*static_cast<MHDFloat>(rand())/RAND_MAX)+this->mMax;
      } else
      {
         srand(this->mStartSeed + seed);
         val = ((this->mMin-this->mMax)*static_cast<MHDFloat>(rand())/RAND_MAX)+this->mMax;
         timespec t;
         clock_gettime(CLOCK_REALTIME, &t);
         srand(t.tv_sec + t.tv_nsec);
      }
   }

   void RandomScalarState::makeRandom(MHDComplex& val, const int i, const int j, const int k) const
   {
      MHDFloat tmp;
      this->makeRandom(tmp, i, j, k);

      val.real() = tmp;

      #if defined QUICC_SPATIALSCHEME_TFT || defined QUICC_SPATIALSCHEME_FFF || defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM || defined QUICC_SPATIALSCHEME_AFT || defined QUICC_SPATIALSCHEME_CFT
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) != 0)
         {
            this->makeRandom(tmp, i, j, k);
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

            this->makeRandom(tmp, i, j, k, seed);
            val.real() = tmp;
            this->makeRandom(tmp, i, j, k, seed + 3);
            if(k2D < n2D/2)
            {
               val.imag() = tmp;
            } else
            {
               val.imag() = -tmp;
            }

         } else if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) != 0)
         {
            this->makeRandom(tmp, i, j, k);
            val.imag() = tmp;
         } else
         {
            val.imag() = 0.0;
         }
      #elif defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k) != 0)
         {
            this->makeRandom(tmp, i, j, k);
            val.imag() = tmp;
         } else
         {
            val.imag() = 0.0;
         }
      #endif //defined QUICC_SPATIALSCHEME_TFT || defined QUICC_SPATIALSCHEME_FFF || defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM || defined QUICC_SPATIALSCHEME_AFT || defined QUICC_SPATIALSCHEME_CFT
   }

}
}
