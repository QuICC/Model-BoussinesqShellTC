/** 
 * @file StreamVerticalWrapper.cpp
 * @brief Source of the streamfunction and vertical velocity to velocity wrapper
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
#include "Diagnostics/StreamVerticalWrapper.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Diagnostics {

   StreamVerticalWrapper::StreamVerticalWrapper(Datatypes::SharedScalarVariableType spStream, Datatypes::SharedScalarVariableType spVertical)
      : mspStream(spStream), mspVertical(spVertical)
   {
   }

   StreamVerticalWrapper::~StreamVerticalWrapper()
   {
   }

   const Datatypes::PhysicalScalarType& StreamVerticalWrapper::one() const
   {
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFF
         // Safety assert
         assert(this->mspVertical);

         return this->mspVertical->dom(0).phys();
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //GEOMHDISCC_SPATIALSCHEME_TFF
   }

   const Datatypes::PhysicalScalarType& StreamVerticalWrapper::two() const
   {
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFF
         // Safety assert
         assert(this->mspStream);

         return this->mspStream->dom(0).grad().comp(FieldComponents::Physical::Y);
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //GEOMHDISCC_SPATIALSCHEME_TFF
   }

   const Datatypes::PhysicalScalarType& StreamVerticalWrapper::three() const
   {
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFF
         // Safety assert
         assert(this->mspStream);

         return this->mspStream->dom(0).grad().comp(FieldComponents::Physical::X);
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //GEOMHDISCC_SPATIALSCHEME_TFF
   }

   const SharedResolution StreamVerticalWrapper::spRes() const
   {
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFF
         // Safety assert
         assert(this->mspStream);

         return this->mspStream->dom(0).spRes();
      #else
         throw Exception("Wrapper is not defined for this scheme");
      #endif //GEOMHDISCC_SPATIALSCHEME_TFF
   }

}
}
