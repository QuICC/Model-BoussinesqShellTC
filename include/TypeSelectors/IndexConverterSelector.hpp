/** \file IndexConverterSelector.hpp
 *  \brief Preprocessor macros to setup the index converters
 */

#ifndef INDEXCONVERTERSELCTOR_HPP
#define INDEXCONVERTERSELCTOR_HPP

//
// Create index converter macro macros
//

// include serial communicators
#include "Communicators/Converters/NoIndexConv.hpp"
#include "Communicators/Converters/SHIndexConv.hpp"

namespace GeoMHDiSCC {

   namespace Parallel {

      template<Dimensions::Transform::Id TId> struct IndexConverterSelector;

      template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
      {
         typedef NoIndexConv  Type;
      };

      template <> struct IndexConverterSelector<Dimensions::Transform::TRA3D>
      {
         typedef NoIndexConv  Type;
      };
   }
}

#endif // INDEXCONVERTERSELCTOR_HPP
