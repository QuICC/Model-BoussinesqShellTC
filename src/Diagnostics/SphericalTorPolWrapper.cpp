/** 
 * @file SphericalTorPolWrapper.cpp
 * @brief Source of the spherical toroidal/poloidal decomposition to velocity wrapper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Debug includes
//
#include <cassert>

// System includes
//

// External includes
//

// Class include
//
#include "Diagnostics/SphericalTorPolWrapper.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Diagnostics {

   SphericalTorPolWrapper::SphericalTorPolWrapper(Datatypes::SharedVectorVariableType spTorPol)
      : mspTorPol(spTorPol)
   {
   }

   SphericalTorPolWrapper::~SphericalTorPolWrapper()
   {
   }

   const Datatypes::PhysicalScalarType& SphericalTorPolWrapper::one() const
   {
      #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFM_TORPOL
         // Safety assert
         assert(this->mspTorPol);

         return this->mspTorPol->dom(0).phys().comp(FieldComponents::Physical::R);
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFM_TORPOL 
   }

   const Datatypes::PhysicalScalarType& SphericalTorPolWrapper::two() const
   {
      #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFM_TORPOL 
         // Safety assert
         assert(this->mspTorPol);

         return this->mspTorPol->dom(0).phys().comp(FieldComponents::Physical::THETA);
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFM_TORPOL 
   }

   const Datatypes::PhysicalScalarType& SphericalTorPolWrapper::three() const
   {
      #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFM_TORPOL 
         // Safety assert
         assert(this->mspTorPol);

         return this->mspTorPol->dom(0).phys().comp(FieldComponents::Physical::PHI);
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFM_TORPOL 
   }

   const SharedResolution SphericalTorPolWrapper::spRes() const
   {
      #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFM_TORPOL 
         // Safety assert
         assert(this->mspTorPol);

         return this->mspTorPol->dom(0).spRes();
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLFM_TORPOL 
   }

}
}
