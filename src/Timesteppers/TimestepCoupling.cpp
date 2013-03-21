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
#include "Exceptions/Exception.hpp"

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
      return std::tr1::get<0>(this->mFields.find(id)->second);
   }

   int TimestepCoupling::idx(const TimestepCoupling::FieldIdType& id) const
   {
      // Safety assert
      assert(this->mFields.count(id) != 0);

      // return the stored index
      return std::tr1::get<1>(this->mFields.find(id)->second);
   }

   void TimestepCoupling::updateType(const int idx, const bool isComplex)
   {
      // Complex flag only changes if isComplex is true
      if(isComplex)
      {
         // Loop over all fields
         std::map<FieldIdType, std::tr1::tuple<bool,int,int> >::iterator it;
         for(it = this->mFields.begin(); it != this->mFields.end(); ++it)
         {
            // Check for the right index
            if(std::tr1::get<1>(it->second) == idx)
            {
               // Set complex flag to true
               std::tr1::get<0>(it->second) = true;
            }
         }
      }
   }

   void TimestepCoupling::checkStart(const int idx, const int start)
   {
      // Loop over all fields
      std::map<FieldIdType, std::tr1::tuple<bool,int,int> >::iterator it;
      for(it = this->mFields.begin(); it != this->mFields.end(); ++it)
      {
         // Check for the right index
         if(std::tr1::get<1>(it->second) == idx)
         {
            if(std::tr1::get<2>(it->second) != start)
            {
               throw Exception("Coupled equations have incompatible starting indexes");
            }
         }
      }
   }

   void TimestepCoupling::updateIndex(const int oldIdx, const int newIdx)
   {
      // Loop over all fields
      std::map<FieldIdType, std::tr1::tuple<bool,int,int> >::iterator it;
      for(it = this->mFields.begin(); it != this->mFields.end(); ++it)
      {
         // Check for the right index
         if(std::tr1::get<1>(it->second) == oldIdx)
         {
            // Set complex flag to true
            std::tr1::get<1>(it->second) = newIdx;
         }
      }
   }

   void TimestepCoupling::addField(const TimestepCoupling::FieldIdType& id, const bool isComplex, const int idx, const int start)
   {
      this->mFields.insert(std::make_pair(id, std::tr1::make_tuple(isComplex, idx, start)));
   }
}
}
