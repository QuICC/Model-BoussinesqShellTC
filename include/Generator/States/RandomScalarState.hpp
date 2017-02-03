/**
 * @file RandomScalarState.hpp
 * @brief Implementation of the equation to generate a random scalar state 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef RANDOMSCALARSTATE_HPP
#define RANDOMSCALARSTATE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate a random scalar state
    */
   class RandomScalarState: public IScalarEquation
   {
      public:
         enum SpecialId {
            ONLYMEAN, ///< Only mean
            NOMEAN, ///< No mean 
            NOTHING, ///< Nothing special
         };

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         RandomScalarState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RandomScalarState();

         /**
          * @brief Compute the random state as a source term
          *
          * @param compId  ID of the spectral component
          * @param i      Fast direction
          * @param j      Middle direction
          * @param k      Slow direction
          */
         virtual Datatypes::SpectralScalarType::PointType sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set spectrum variables
          */
         void setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const SpecialId id = NOTHING);

         /**
          * @brief Set spectrum variables
          */
         void setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D, const SpecialId id = NOTHING);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

         /**
          * Generate Random value
          */
         void makeRandom(MHDFloat& val, const int i1D, const int i3D, const int i2D, const unsigned int seed = 1) const;

         /**
          * Generate Random value
          */
         void makeRandom(MHDComplex& val, const int i1D, const int i3D, const int i2D) const;

      private:
         /**
          * @brief Minimum value in coefficient range
          */
         MHDFloat mMin;

         /**
          * @brief Maximum value in coefficient range
          */
         MHDFloat mMax;

         /**
          * @brief Ratio between first and last coefficient in first direction
          */
         MHDFloat mRatio1D;

         /**
          * @brief Ratio between first and last coefficient in second direction
          */
         MHDFloat mRatio2D;

         /**
          * @brief Ratio between first and last coefficient in third direction
          */
         MHDFloat mRatio3D;

         /**
          * @brief Starting seed used for random generator
          */
         int mStartSeed;

         /**
          * @brief ID flag for special setup
          */
         SpecialId mSpecial;
   };

   /// Typedef for a shared RandomScalarState
   typedef SharedPtrMacro<RandomScalarState> SharedRandomScalarState;

}
}

#endif // RANDOMSCALARSTATE_HPP
