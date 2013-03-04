/** \file TimestepCoupling.cpp
 *  \brief Source of the implementation of the equation coupling at the timestepper level
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/TimestepCoupling.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Timestep {

   TimestepCoupling::TimestepCoupling()
      : mCounter(0)
   {
   }

   TimestepCoupling::~TimestepCoupling()
   {
   }

   bool TimestepCoupling::isPresent(const TimestepCoupling::FieldIdType& id) const
   {
      return this->mFields.count(id);
   }

   int TimestepCoupling::newIndex()
   {
      // Increment the counter
      ++this->mCounter;

      // return a negative index
      return -this->mCounter;
   }

   bool TimestepCoupling::isComplex(const TimestepCoupling::FieldIdType& id) const
   {
      // Safety assert
      assert(this->mFields.count(id) != 0);

      // return the stored index
      return this->mFields.find(id)->second.first;
   }

   int TimestepCoupling::idx(const TimestepCoupling::FieldIdType& id) const
   {
      // Safety assert
      assert(this->mFields.count(id) != 0);

      // return the stored index
      return this->mFields.find(id)->second.second;
   }

   void TimestepCoupling::updateType(int idx, bool isComplex)
   {
      // Complex flag only changes if isComplex is true
      if(isComplex)
      {
         // Loop over all fields
         std::map<FieldIdType, std::pair<bool,int> >::iterator it;
         for(it = this->mFields.begin(); it != this->mFields.end(); ++it)
         {
            // Check for the right index
            if(it->second.second == idx)
            {
               // Set complex flag to true
               it->second.first = true;
            }
         }
      }
   }

   void TimestepCoupling::updateIndex(int oldIdx, int newIdx)
   {
      // Loop over all fields
      std::map<FieldIdType, std::pair<bool,int> >::iterator it;
      for(it = this->mFields.begin(); it != this->mFields.end(); ++it)
      {
         // Check for the right index
         if(it->second.second == oldIdx)
         {
            // Set complex flag to true
            it->second.second = newIdx;
         }
      }
   }

   void TimestepCoupling::addField(const TimestepCoupling::FieldIdType& id, bool isComplex, int idx)
   {
      this->mFields.insert(std::make_pair(id, std::make_pair(isComplex, idx)));
   }
}
}
