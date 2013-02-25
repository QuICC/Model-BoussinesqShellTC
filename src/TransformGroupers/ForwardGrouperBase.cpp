/** \file ForwardGrouperBase.cpp
 *  \brief Source of the implementation of the equation wise forward transform grouper
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
#include "TransformGroupers/ForwardGrouperBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   ForwardGrouperBase::ForwardGrouperBase()
      : split(Splittings::Locations::NONE), mScalarPacks1D(1), mScalarPacks2D(1), mVectorPacks1D(3), mVectorPacks2D(3)
   {
   }

   ForwardGrouperBase::~ForwardGrouperBase()
   {
   }

   ArrayI ForwardGrouperBase::listPacks1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   {
      // Create list of packet sizes
      std::set<int>  list;

      // loop over all variable information
      std::map<PhysicalNames::Id, std::pair<bool,TriBool> >::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange
         if(infoIt->second.second(0))
         {
            if(infoIt->second.first)
            {
               list.insert(this->mScalarPacks1D);
            } else
            {
               list.insert(this->mVectorPacks1D);
            }
         }
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

   ArrayI ForwardGrouperBase::listPacks2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   {  
      // Create list of packet sizes
      std::set<int>  list;

      // loop over all variable information
      std::map<PhysicalNames::Id, std::pair<bool,TriBool> >::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physica field packs for second exchange
         if(infoIt->second.second(0))
         {
            if(infoIt->second.first)
            {
               list.insert(this->mScalarPacks2D);
            } else
            {
               list.insert(this->mVectorPacks2D);
            }
         }
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

   ArrayI ForwardGrouperBase::groupPacks1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   {
      // Initialise the number of packs
      ArrayI packs(1);
      packs.setConstant(0);

      // loop over all variable information
      std::map<PhysicalNames::Id, std::pair<bool,TriBool> >::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange
         if(infoIt->second.second(0))
         {
            if(infoIt->second.first)
            {
               packs(0) += this->mScalarPacks1D;
            } else
            {
               packs(0) += this->mVectorPacks1D;
            }
         }
      }

      return packs;
   }

   ArrayI ForwardGrouperBase::groupPacks2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   {  
      // Initialise the number of packs
      ArrayI packs(1);
      packs.setConstant(0);

      // loop over all variable information
      std::map<PhysicalNames::Id, std::pair<bool,TriBool> >::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for second exchange
         if(infoIt->second.second(0))
         {
            if(infoIt->second.first)
            {
               packs(0) += this->mScalarPacks2D;
            } else
            {
               packs(0) += this->mVectorPacks2D;
            }
         }
      }

      return packs;
   }

}
}
