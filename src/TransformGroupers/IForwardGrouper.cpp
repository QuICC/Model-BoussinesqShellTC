/** 
 * @file IForwardGrouper.cpp
 * @brief Source of the implementation of the equation wise forward transform grouper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <set>

// Class include
//
#include "TransformGroupers/IForwardGrouper.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   IForwardGrouper::IForwardGrouper()
      : split(Splitting::Locations::NONE), mcScalarPacks1D(1), mcScalarPacks2D(1), mcVectorPacks1D(3), mcVectorPacks2D(3)
   {
   }

   IForwardGrouper::~IForwardGrouper()
   {
   }

   ArrayI IForwardGrouper::namePacks1D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {
      // Create list of packet sizes
      std::set<int>  list;

      int count = 0;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange
         if(infoIt->second.needSpectral() && nonInfo.count(infoIt->first) > 0)
         {
            if(infoIt->second.isScalar())
            {
               count = this->mcScalarPacks1D;
            } else
            {
               count = this->mcVectorPacks1D;
            }

            list.insert(count);
         }

         // Add named pack size
         this->mNamedPacks1D.insert(std::make_pair(infoIt->first, count));

         // reset counter
         count = 0;
      }

      // Initialise the number of packs
      ArrayI packs(list.size());

      // Set packet sizes
      std::set<int>::iterator it;
      int i = 0;
      for(it = list.begin(); it != list.end(); ++it, ++i)
      {
         packs(i) = *it;
      }

      return packs;
   }

   ArrayI IForwardGrouper::namePacks2D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {  
      // Create list of packet sizes
      std::set<int>  list;

      int count = 0;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for second exchange
         if(infoIt->second.needSpectral() && nonInfo.count(infoIt->first) > 0)
         {
            if(infoIt->second.isScalar())
            {
               count = this->mcScalarPacks2D;
            } else
            {
               count = this->mcVectorPacks2D;
            }

            list.insert(count);
         }

         // Add named pack size
         this->mNamedPacks2D.insert(std::make_pair(infoIt->first, count));

         // reset counter
         count = 0;
      }

      // Initialise the number of packs
      ArrayI packs(list.size());

      // Set packet sizes
      std::set<int>::iterator it;
      int i = 0;
      for(it = list.begin(); it != list.end(); ++it, ++i)
      {
         packs(i) = *it;
      }

      return packs;
   }

   ArrayI IForwardGrouper::groupPacks1D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {
      // Initialise the number of packs
      ArrayI packs(1);
      packs.setConstant(0);

      // Compute some of packs
      packs(0) = this->namePacks1D(varInfo, nonInfo).sum();

      return packs;
   }

   ArrayI IForwardGrouper::groupPacks2D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {  
      // Initialise the number of packs
      ArrayI packs(1);
      packs.setConstant(0);

      // Compute some of packs
      packs(0) = this->namePacks2D(varInfo, nonInfo).sum();

      return packs;
   }

}
}
