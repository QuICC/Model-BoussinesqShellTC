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
      : split(Splitting::Locations::NONE), mcScalarPacks1D(1), mcScalarPacks2D(1), mcVectorPacks1D(3), mcVectorPacks2D(3), mScalarPacks1D(0), mScalarPacks2D(0), mVectorPacks1D(0), mVectorPacks2D(0)
   {
   }

   IForwardGrouper::~IForwardGrouper()
   {
   }

   ArrayI IForwardGrouper::listPacks1D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {
      // Create list of packet sizes
      std::set<int>  list;

      int countScalar = 0;
      int countVector = 0;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange
         if(infoIt->second.needSpectral() && nonInfo.count(infoIt->first) > 0)
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mcScalarPacks1D);
               ++countScalar;
            } else
            {
               list.insert(this->mcVectorPacks1D);
               ++countVector;
            }
         }
      }

      // Set to zero if no transform is involved
      if(countScalar > 0)
      {
         this->mScalarPacks1D = this->mcScalarPacks1D;
      } else
      {
         this->mScalarPacks1D = 0;
      }

      // Set to zero if no transform is involved
      if(countVector > 0)
      {
         this->mVectorPacks1D = this->mcVectorPacks1D;
      } else
      {
         this->mVectorPacks1D = 0;
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

   ArrayI IForwardGrouper::listPacks2D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {  
      // Create list of packet sizes
      std::set<int>  list;

      int countScalar = 0;
      int countVector = 0;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for second exchange
         if(infoIt->second.needSpectral() && nonInfo.count(infoIt->first) > 0)
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mcScalarPacks2D);
               ++countScalar;
            } else
            {
               list.insert(this->mcVectorPacks2D);
               ++countVector;
            }
         }
      }

      // Set to zero if no transform is involved
      if(countScalar > 0)
      {
         this->mScalarPacks2D = this->mcScalarPacks2D;
      } else
      {
         this->mScalarPacks2D = 0;
      }

      // Set to zero if no transform is involved
      if(countVector > 0)
      {
         this->mVectorPacks2D = this->mcVectorPacks2D;
      } else
      {
         this->mVectorPacks2D = 0;
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

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical -> spectral field packs for first exchange
         if(infoIt->second.needSpectral() && nonInfo.count(infoIt->first) > 0)
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

   ArrayI IForwardGrouper::groupPacks2D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {  
      // Initialise the number of packs
      ArrayI packs(1);
      packs.setConstant(0);

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical -> spectral field packs for second exchange
         if(infoIt->second.needSpectral() && nonInfo.count(infoIt->first) > 0)
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
