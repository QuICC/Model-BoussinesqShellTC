/** 
 * @file ISchemeCosts.cpp
 * @brief Source of the base for a cost based scheme implementations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "SpatialSchemes/ISchemeCosts.hpp"

// Project includes
//

namespace QuICC {

namespace Schemes {

   ISchemeCosts::ISchemeCosts(const int dims)
      : mCosts(dims), mScalings(dims), mMemory(dims)
   {
   }

   ISchemeCosts::~ISchemeCosts()
   {
   }

   void ISchemeCosts::resetLoad()
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

   void ISchemeCosts::updateLoad(const int parts)
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

   Array ISchemeCosts::loadWeights()
   {
      return this->mCosts.array() * this->mScalings.array();
   }

   double ISchemeCosts::memoryScore(SharedResolution spRes)
   {
      #ifdef QUICC_MEMORYUSAGE_LIMITED
         return this->mMemory.prod();
      #else
         return 1.0;
      #endif //QUICC_MEMORYUSAGE_LIMITED
   }

   void ISchemeCosts::setCost(MHDFloat c, Dimensions::Transform::Id id)
   {
      // Assert for positive cost
      assert(c > 0);

      // Assert for number of transforms
      assert(static_cast<int>(id) < this->mCosts.size());

      this->mCosts(static_cast<int>(id)) = c;
   }

   void ISchemeCosts::setScaling(MHDFloat c, Dimensions::Transform::Id id)
   {
      // Assert for positive cost
      assert(c > 0);

      // Assert for number of transforms
      assert(static_cast<int>(id) < this->mScalings.size());

      this->mScalings(static_cast<int>(id)) = c;
   }

   void ISchemeCosts::setMemory(MHDFloat c, Dimensions::Transform::Id id)
   {
      // Assert for positive cost
      assert(c > 0);

      // Assert for number of transforms
      assert(static_cast<int>(id) < this->mMemory.size());

      this->mMemory(static_cast<int>(id)) = c;
   }

}
}
