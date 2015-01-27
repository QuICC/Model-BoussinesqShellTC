/** 
 * @file IndexCounter.cpp
 * @brief Source of base class for the index counters
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Resolutions/Tools/IndexCounter.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   IndexCounter::IndexCounter()
   {
   }

   IndexCounter::~IndexCounter()
   {
   }

   std::tr1::tuple<int,int,int> IndexCounter::makeKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const
   {
      std::tr1::tuple<int,int,int> key;

      if(id == Dimensions::Transform::TRA1D)
      {
         key = std::tr1::make_tuple(i, k, j);
      } else if(id == Dimensions::Transform::TRA2D)
      {
         key = std::tr1::make_tuple(j, i, k);
      } else if(id == Dimensions::Transform::TRA3D)
      {
         key = std::tr1::make_tuple(k, j, i);
      }

      return key;
   }

}
