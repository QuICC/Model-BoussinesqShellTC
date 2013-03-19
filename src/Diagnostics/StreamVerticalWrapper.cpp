/** \file StreamVerticalWrapper.cpp
 *  \brief Source of the streamfunction and vertical velocity to velocity wrapper
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
      // Safety assert
      assert(this->mspStream);

      return this->mspStream->dom(0).grad().comp(FieldComponents::Physical::TWO);
   }

   const Datatypes::PhysicalScalarType& StreamVerticalWrapper::two() const
   {
      // Safety assert
      assert(this->mspStream);

      return this->mspStream->dom(0).grad().comp(FieldComponents::Physical::ONE);
   }

   const Datatypes::PhysicalScalarType& StreamVerticalWrapper::three() const
   {
      // Safety assert
      assert(this->mspVertical);

      return this->mspVertical->dom(0).phys();
   }

}
}
