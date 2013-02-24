/** \file SchemeBase.cpp
 *  \brief Source of the base for a cost based scheme implementations
 */

// System includes
//

// External includes
//

// Class include
//
#include "Base/SpatialSchemes/SchemeBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SchemeBase::SchemeBase(const int dims)
      : mCosts(dims), mScalings(dims), mMemory(dims)
   {
   }

   SchemeBase::~SchemeBase()
   {
   }

   void SchemeBase::resetLoad()
   {
      // Clear load list
      this->mLoadList.clear();

      // Clear optimal loads
      while (!this->mOptimalLoad.empty())
      {
         this->mOptimalLoad.pop();
      }

      // Clear load sums
      this->mLoad.clear();
      
      // Clear regularized loads per part
      this->mRegularLoad.clear();
   }

   void SchemeBase::updateLoad(const int parts)
   {
      // Create typedef to simplify notation
      typedef  std::multimap<int, int>::const_iterator MapIt;

      // Const iterator to map object
      MapIt it;

      // Pair of iterator to work as range
      std::pair<MapIt, MapIt> range;

      // Loop over all parts
      for(int i = 0; i < parts; ++i)
      {
         // Reset loads
         this->mLoad.at(i) = 0;

         // Get assigned loads
         range = this->mRegularLoad.equal_range(i);

         // Loop over the assigned loads
         for(it = range.first; it != range.second; ++it)
         {
            this->mLoad.at(i) += it->second;
         }
      }
   }

}
