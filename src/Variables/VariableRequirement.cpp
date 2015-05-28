/** 
 * @file VariableRequirement.cpp
 * @brief Source of the variable requirements
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "Variables/VariableRequirement.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

   VariableRequirement::VariableRequirement()
      : mNoField(false,false,false,false)
   {
   }

   VariableRequirement::~VariableRequirement()
   {
   }

   const FieldRequirement& VariableRequirement::field(const PhysicalNames::Id id) const
   {
      if(this->mInfo.count(id) == 0)
      {
         return this->mNoField;
      } else
      {
         return this->mInfo.find(id)->second;
      }
   }

   FieldRequirement& VariableRequirement::rField(const PhysicalNames::Id id)
   {
      if(this->mInfo.count(id) == 0)
      {
         throw Exception("Tried to modify requirements of inexistant field!");
      }

      return this->mInfo.find(id)->second;
   }

   void VariableRequirement::addField(const PhysicalNames::Id id, const FieldRequirement& req)
   {
      this->mInfo.insert(std::make_pair(id,req));
   }

   VariableRequirement::const_iterator VariableRequirement::begin() const
   {
      return this->mInfo.begin();
   }

   VariableRequirement::const_iterator VariableRequirement::end() const
   {
      return this->mInfo.end();
   }

   void VariableRequirement::merge(const VariableRequirement& req)
   {
      // Create iterator
      VariableRequirement::const_iterator it;

      for(it = req.begin(); it != req.end(); it++)
      {
         if(this->mInfo.count(it->first) == 0)
         {
            this->mInfo.insert(*it);
         } else
         {
            this->mInfo.find(it->first)->second.merge(it->second);
         }
      }
   }
}
