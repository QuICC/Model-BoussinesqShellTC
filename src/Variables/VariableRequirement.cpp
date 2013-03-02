/** \file VariableRequirement.cpp
 *  \brief Source of the variable requirements
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

namespace GeoMHDiSCC {

   VariableRequirement::VariableRequirement()
   {
   }

   VariableRequirement::~VariableRequirement()
   {
   }

   const FieldRequirement& VariableRequirement::requirement(const PhysicalNames::Id id) const
   {
      // Safety assert
      assert(this->mInfo.find(id) != this->mInfo.end());

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
}
