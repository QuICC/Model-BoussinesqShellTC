/** 
 * @file CartesianTorPolWrapper.cpp
 * @brief Source of the cartesian toroidal/poloidal decomposition to velocity wrapper
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
#include "Diagnostics/CartesianTorPolWrapper.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Diagnostics {

   CartesianTorPolWrapper::CartesianTorPolWrapper(Datatypes::SharedVectorVariableType spTorPol)
      : mspTorPol(spTorPol)
   {
   }

   CartesianTorPolWrapper::~CartesianTorPolWrapper()
   {
   }

   const Datatypes::PhysicalScalarType& CartesianTorPolWrapper::one() const
   {
      #if defined GEOMHDISCC_SPATIALSCHEME_TFF_TORPOL
         // Safety assert
         assert(this->mspTorPol);

         return this->mspTorPol->dom(0).phys().comp(FieldComponents::Physical::Z);
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //defined GEOMHDISCC_SPATIALSCHEME_TFF_TORPOL 
   }

   const Datatypes::PhysicalScalarType& CartesianTorPolWrapper::two() const
   {
      #if defined GEOMHDISCC_SPATIALSCHEME_TFF_TORPOL 
         // Safety assert
         assert(this->mspTorPol);

         return this->mspTorPol->dom(0).phys().comp(FieldComponents::Physical::X);
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //defined GEOMHDISCC_SPATIALSCHEME_TFF_TORPOL
   }

   const Datatypes::PhysicalScalarType& CartesianTorPolWrapper::three() const
   {
      #if defined GEOMHDISCC_SPATIALSCHEME_TFF_TORPOL
         // Safety assert
         assert(this->mspTorPol);

         return this->mspTorPol->dom(0).phys().comp(FieldComponents::Physical::Y);
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //defined GEOMHDISCC_SPATIALSCHEME_TFF_TORPOL
   }

   const SharedResolution CartesianTorPolWrapper::spRes() const
   {
      #if defined GEOMHDISCC_SPATIALSCHEME_TFF_TORPOL
         // Safety assert
         assert(this->mspTorPol);

         return this->mspTorPol->dom(0).spRes();
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //defined GEOMHDISCC_SPATIALSCHEME_TFF_TORPOL
   }

}
}
