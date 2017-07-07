/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq spherical Couette around
 * an equatorially sprecessing orthogonal axis in a spherical shell model
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//
#include "TypeSelectors/TransformSelector.hpp"

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Shell/NuttatingCouette/Momentum.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/Cross.hpp"
#include "PhysicalOperators/SphericalCoriolis.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace NuttatingCouette {

   Momentum::Momentum(SharedEquationParameters spEqParams)
      : Couette::MomentumBase(spEqParams)
   {
   }

   Momentum::~Momentum()
   {
   }

   Datatypes::SpectralScalarType::PointType Momentum::boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {  
      if(compId == FieldComponents::Spectral::TOR)
      {
         #if defined QUICC_SPATIALSCHEME_SLFL
            if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 1)
            {
               if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 1)
               {
                  if(i == 1)
                  {
         #elif defined QUICC_SPATIALSCHEME_SLFM
            if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 1)
            {
               if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 1)
               {
                  if(i == 1)
                  {
         #endif //defined QUICC_SPATIALSCHEME_SLFL

                     MHDFloat sgn = std::pow(-1.0,std::signbit(this->eqParams().nd(NonDimensional::ROSSBY)));
                     MHDFloat ri = this->eqParams().nd(NonDimensional::RO)*this->eqParams().nd(NonDimensional::RRATIO);
                     MHDFloat norm = std::sqrt(3.0/(4.0*Math::PI));
                     MHDFloat omega = this->eqParams().nd(NonDimensional::OMEGA);
					 MHDFloat t = this->time();
                     MHDComplex phase = std::exp(Math::cI * omega * t);
                     return Datatypes::SpectralScalarType::PointType(phase*sgn*ri/norm);
                  }
               }
            }

         return Datatypes::SpectralScalarType::PointType(0.0);
      } else
      {
         throw std::logic_error("Poloidal component should not set inhomogeneous boundary condition");

         return Datatypes::SpectralScalarType::PointType(0.0);
      }
   }

}
}
}
}
}
