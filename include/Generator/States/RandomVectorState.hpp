/**
 * @file RandomVectorState.hpp
 * @brief Implementation of the equation to generate a random vector state 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef RANDOMVECTORSTATE_HPP
#define RANDOMVECTORSTATE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate a random vector state
    */
   class RandomVectorState: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         RandomVectorState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RandomVectorState();

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
         void setSpectrum(const FieldComponents::Spectral::Id compId, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D);

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
          * @brief Set the nonliner integration components
          */
         virtual void setNLComponents();

         /**
          * Generate Random value
          */
         void makeRandom(MHDFloat& val, const int i1D, const int i3D, const int i2D, const MHDFloat minVal, const MHDFloat maxVal) const;

         /**
          * Generate Random value
          */
         void makeRandom(MHDComplex& val, const int i1D, const int i3D, const int i2D, const MHDFloat minVal, const MHDFloat maxVal) const;

      private:
         /**
          * @brief Minimum value in coefficient range
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mMin;

         /**
          * @brief Maximum value in coefficient range
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mMax;

         /**
          * @brief Ratio between first and last coefficient in first direction
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mRatio1D;

         /**
          * @brief Ratio between first and last coefficient in second direction
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mRatio2D;

         /**
          * @brief Ratio between first and last coefficient in third direction
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mRatio3D;
   };

   /// Typedef for a shared RandomVectorState
   typedef SharedPtrMacro<RandomVectorState> SharedRandomVectorState;

}
}

#endif // RANDOMVECTORSTATE_HPP
