/** \file BackwardGrouperBase.cpp
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
#include "TransformGroupers/BackwardGrouperBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   BackwardGrouperBase::BackwardGrouperBase()
      : split(Splitting::Locations::NONE), mScalarPacks1D(1), mGradientPacks1D(2), mVectorPacks1D(3), mCurlPacks1D(3), mScalarPacks2D(1), mGradientPacks2D(2), mVectorPacks2D(3), mCurlPacks2D(3)
   {
   }

   BackwardGrouperBase::~BackwardGrouperBase()
   {
   }

   ArrayI BackwardGrouperBase::listPacks1D(const VariableRequirement& varInfo)
   {
      // Create list of packet sizes
      std::set<int>  list;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mScalarPacks1D);
            } else
            {
               list.insert(this->mVectorPacks1D);
            }
         }

         // add physical differential field packs first exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mGradientPacks1D);
            } else
            {
               list.insert(this->mCurlPacks1D);
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

   ArrayI BackwardGrouperBase::listPacks2D(const VariableRequirement& varInfo)
   {  
      // Create list of packet sizes
      std::set<int>  list;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for second exchange
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mScalarPacks2D);
            } else
            {
               list.insert(this->mVectorPacks2D);
            }
         }

         // add physical differential field packs for second exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mGradientPacks2D);
            } else
            {
               list.insert(this->mCurlPacks2D);
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

   ArrayI BackwardGrouperBase::namePacks1D(const VariableRequirement& varInfo)
   {
      // Create list of packet sizes
      std::set<int>  list;
      int tmp = 0;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               tmp += this->mScalarPacks1D;
            } else
            {
               tmp += this->mVectorPacks1D;
            }
         }

         // add physical differential field packs for first exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               tmp += this->mGradientPacks1D;
            } else
            {
               tmp += this->mCurlPacks1D;
            }
         }

         // Add in the combine pack sizes
         if(tmp > 0)
         {
            // Add combined pack size to list
            list.insert(tmp);

            // Add named pack size
            this->mNamedPacks1D.insert(std::make_pair(infoIt->first, tmp));

            // reset counter
            tmp = 0;
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

   ArrayI BackwardGrouperBase::namePacks2D(const VariableRequirement& varInfo)
   {  
      // Create list of packet sizes
      std::set<int>  list;
      int tmp = 0;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for second exchange
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               tmp += this->mScalarPacks2D;
            } else
            {
               tmp += this->mVectorPacks2D;
            }
         }

         // add physical differential field packs for second exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               tmp += this->mGradientPacks2D + 1;
            } else
            {
               tmp += this->mCurlPacks2D;
            }
         }

         // Add in the combine pack sizes
         if(tmp > 0)
         {
            // Add combined pack size to list
            list.insert(tmp);

            // Add named pack size
            this->mNamedPacks2D.insert(std::make_pair(infoIt->first, tmp));

            // reset counter
            tmp = 0;
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

   ArrayI BackwardGrouperBase::groupPacks1D(const VariableRequirement& varInfo)
   {
      // Initialise the number of packs
      ArrayI packs(1);
      packs.setConstant(0);

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange 
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               packs(0) += this->mScalarPacks1D;
            } else
            {
               packs(0) += this->mVectorPacks1D;
            }
         }

         // add physical differential field packs for first exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               packs(0) += this->mGradientPacks1D;
            } else
            {
               packs(0) += this->mCurlPacks1D;
            }
         }
      }

      return packs;
   }

   ArrayI BackwardGrouperBase::groupPacks2D(const VariableRequirement& varInfo)
   {  
      // Initialise the number of packs
      ArrayI packs(1);
      packs.setConstant(0);

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for second exchange
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               packs(0) += this->mScalarPacks2D;
            } else
            {
               packs(0) += this->mVectorPacks2D;
            }
         }

         // add physical differential field packs for second exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               packs(0) += this->mGradientPacks2D + 1;
            } else
            {
               packs(0) += this->mCurlPacks2D;
            }
         }
      }

      return packs;
   }

}
}
