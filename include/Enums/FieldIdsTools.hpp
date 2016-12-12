/**
 * @file FieldIds.hpp
 * @brief Definition of some useful function to work with field IDs 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FIELDIDSTOOLS_HPP
#define FIELDIDSTOOLS_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <utility>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"

namespace GeoMHDiSCC {
   
   /**
    * Convert field component to array index
    */
   int fieldIndex(const FieldComponents::Physical::Id id);
   
   /**
    * Convert field component pair to symmetric matrix indexes pair
    */
   std::pair<int,int> fieldPairSym(const FieldComponents::Physical::Id id1, const FieldComponents::Physical::Id id2);
}

#endif // FIELDIDSTOOLS_HPP
