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

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate a random scalar state
    */
   class RandomScalarState: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param pyName     Python script name
          * @param spEqParams Shared equation parameters
          */
         RandomScalarState(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RandomScalarState();

         /**
          * @brief Compute the random state as a source term
          *
          * @param compId  ID of the spectral component
          * @param iX      Index for the X direction
          * @param iZ      Index for the Z direction
          * @param iY      Index for the Y direction
          */
         virtual Datatypes::SpectralScalarType::PointType sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set spectrum variables
          */
         void setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat xRatio, const MHDFloat yRatio, const MHDFloat zRatio);

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
         void makeRandom(MHDFloat& val, const int iX, const int iZ, const int iY) const;

         /**
          * Generate Random value
          */
         void makeRandom(MHDComplex& val, const int iX, const int iZ, const int iY) const;

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
          * @brief Ratio between first and last coefficient in X direction
          */
         MHDFloat mXRatio;

         /**
          * @brief Ratio between first and last coefficient in Y direction
          */
         MHDFloat mYRatio;

         /**
          * @brief Ratio between first and last coefficient in Z direction
          */
         MHDFloat mZRatio;
   };

   /// Typedef for a shared RandomScalarState
   typedef SharedPtrMacro<RandomScalarState> SharedRandomScalarState;

}
}

#endif // RANDOMSCALARSTATE_HPP
