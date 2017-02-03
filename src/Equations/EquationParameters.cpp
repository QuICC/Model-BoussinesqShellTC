/** 
 * @file EquationParameters.cpp
 * @brief Source of the implementation of the non dimensional parameters
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/EquationParameters.hpp"

// Project includes
//
#include "IoTools/IdToHuman.hpp"
#include "IoTools/HumanToId.hpp"

namespace QuICC {

namespace Equations {

   EquationParameters::EquationParameters()
   {
   }

   EquationParameters::~EquationParameters()
   {
   }

   MHDFloat EquationParameters::nd(NonDimensional::Id name) const
   {
      return this->mND.at(name);
   }

   const std::map<NonDimensional::Id,MHDFloat>& EquationParameters::map() const
   {
      return this->mND;
   }

   std::map<NonDimensional::Id,MHDFloat>& EquationParameters::rMap()
   {
      return this->mND;
   }

   std::vector<NonDimensional::Id> EquationParameters::ids() const
   {
      // Storage for the IDs
      std::vector<NonDimensional::Id> ids;

      std::map<NonDimensional::Id,MHDFloat>::const_iterator it;
      for(it = this->mND.begin(); it != this->mND.end(); it++)
      {
         ids.push_back(it->first);
      }

      return ids;
   }

   std::vector<std::string> EquationParameters::names() const
   {
      // Storage for the names
      std::vector<std::string> names;

      std::map<NonDimensional::Id,MHDFloat>::const_iterator it;
      for(it = this->mND.begin(); it != this->mND.end(); it++)
      {
         names.push_back(IoTools::IdToHuman::toTag(it->first));
      }

      return names;
   }

   void EquationParameters::init(const std::map<std::string,MHDFloat>& parameters)
   {
      // Add non dimensional parameters
      std::map<std::string,MHDFloat>::const_iterator it;
      for(it = parameters.begin(); it != parameters.end(); it++)
      {
         this->mND.insert(std::make_pair(IoTools::HumanToId::toNd(it->first), it->second));
      }
   }

}
}
