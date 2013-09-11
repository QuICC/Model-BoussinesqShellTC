/** 
 * @file IBackwardGrouper.cpp
 * @brief Source of the implementation of the equation wise forward transform grouper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
#include "TransformGroupers/IBackwardGrouper.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   IBackwardGrouper::IBackwardGrouper()
      : split(Splitting::Locations::NONE), mcScalarPacks1D(1), mcGradientPacks1D(2), mcVectorPacks1D(3), mcCurlPacks1D(3), mcScalarPacks2D(1), mcGradientPacks2D(2), mcVectorPacks2D(3), mcCurlPacks2D(3)
   {
   }

   IBackwardGrouper::~IBackwardGrouper()
   {
   }

   ArrayI IBackwardGrouper::listPacks1D(const VariableRequirement& varInfo)
   {
      // Create list of packet sizes
      std::set<int>  list;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add spectral -> physical field packs for first exchange
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mcScalarPacks1D);
            } else
            {
               list.insert(this->mcVectorPacks1D);
            }
         }

         // add spectral -> physical differential field packs first exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mcGradientPacks1D);
            } else
            {
               list.insert(this->mcCurlPacks1D);
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

   ArrayI IBackwardGrouper::listPacks2D(const VariableRequirement& varInfo)
   {  
      // Create list of packet sizes
      std::set<int>  list;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add spectral -> physical field packs for second exchange
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mcScalarPacks2D);
            } else
            {
               list.insert(this->mcVectorPacks2D);
            }
         }

         // add spectral -> physical differential field packs for second exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               list.insert(this->mcGradientPacks2D);
            } else
            {
               list.insert(this->mcCurlPacks2D);
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

   ArrayI IBackwardGrouper::namePacks1D(const VariableRequirement& varInfo)
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
               tmp += this->mcScalarPacks1D;
            } else
            {
               tmp += this->mcVectorPacks1D;
            }
         }

         // add physical differential field packs for first exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               tmp += this->mcGradientPacks1D;
            } else
            {
               tmp += this->mcCurlPacks1D;
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

   ArrayI IBackwardGrouper::namePacks2D(const VariableRequirement& varInfo)
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
               tmp += this->mcScalarPacks2D;
            } else
            {
               tmp += this->mcVectorPacks2D;
            }
         }

         // add physical differential field packs for second exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               tmp += this->mcGradientPacks2D + 1;
            } else
            {
               tmp += this->mcCurlPacks2D;
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

   ArrayI IBackwardGrouper::groupPacks1D(const VariableRequirement& varInfo)
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
               packs(0) += this->mcScalarPacks1D;
            } else
            {
               packs(0) += this->mcVectorPacks1D;
            }
         }

         // add physical differential field packs for first exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               packs(0) += this->mcGradientPacks1D;
            } else
            {
               packs(0) += this->mcCurlPacks1D;
            }
         }
      }

      return packs;
   }

   ArrayI IBackwardGrouper::groupPacks2D(const VariableRequirement& varInfo)
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
               packs(0) += this->mcScalarPacks2D;
            } else
            {
               packs(0) += this->mcVectorPacks2D;
            }
         }

         // add physical differential field packs for second exchange
         if(infoIt->second.needPhysicalDiff())
         {
            if(infoIt->second.isScalar())
            {
               packs(0) += this->mcGradientPacks2D + 1;
            } else
            {
               packs(0) += this->mcCurlPacks2D;
            }
         }
      }

      return packs;
   }

}
}
