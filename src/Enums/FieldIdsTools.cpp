/** 
 * @file FieldIdsTools.cpp
 * @brief Source of useful functions to work with field IDs
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Enums/FieldIdsTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   int fieldIndex(FieldComponents::Physical::Id id)
   {
      if(id == FieldComponents::Physical::ONE)
      {
         return 0;

      } else if(id == FieldComponents::Physical::TWO) 
      {
         return 1;

      } else 
      {
         return 2;
      }
   }

   std::pair<int,int> fieldPairSym(FieldComponents::Physical::Id id1, FieldComponents::Physical::Id id2)
   {
      if(id1 == FieldComponents::Physical::ONE)
      {
         return std::make_pair(fieldIndex(id1), fieldIndex(id2));

      } else if(id1 == FieldComponents::Physical::TWO && id2 != FieldComponents::Physical::ONE) 
      {
         return std::make_pair(fieldIndex(id1), fieldIndex(id2));

      } else if(id1 == FieldComponents::Physical::THREE && id1 == id2) 
      {
         return std::make_pair(fieldIndex(id1), fieldIndex(id2));
      } else
      {
         return std::make_pair(fieldIndex(id2), fieldIndex(id1));
      }
   }

}
