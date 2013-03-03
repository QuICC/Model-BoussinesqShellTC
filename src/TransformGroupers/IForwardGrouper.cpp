/** \file IForwardGrouper.cpp
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
#include "TransformGroupers/IForwardGrouper.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   IForwardGrouper::IForwardGrouper()
      : split(Splitting::Locations::NONE), mScalarPacks1D(1), mScalarPacks2D(1), mVectorPacks1D(3), mVectorPacks2D(3)
   {
   }

   IForwardGrouper::~IForwardGrouper()
   {
   }

   ArrayI IForwardGrouper::listPacks1D(const VariableRequirement& varInfo)
   {
      // Create list of packet sizes
      std::set<int>  list;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange
         if(infoIt->second.needSpectral())
         {
            if(infoIt->second.isScalar())
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

   ArrayI IForwardGrouper::listPacks2D(const VariableRequirement& varInfo)
   {  
      // Create list of packet sizes
      std::set<int>  list;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physica field packs for second exchange
         if(infoIt->second.needSpectral())
         {
            if(infoIt->second.isScalar())
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

   ArrayI IForwardGrouper::groupPacks1D(const VariableRequirement& varInfo)
   {
      // Initialise the number of packs
      ArrayI packs(1);
      packs.setConstant(0);

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange
         if(infoIt->second.needSpectral())
         {
            if(infoIt->second.isScalar())
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

   ArrayI IForwardGrouper::groupPacks2D(const VariableRequirement& varInfo)
   {  
      // Initialise the number of packs
      ArrayI packs(1);
      packs.setConstant(0);

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for second exchange
         if(infoIt->second.needSpectral())
         {
            if(infoIt->second.isScalar())
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
