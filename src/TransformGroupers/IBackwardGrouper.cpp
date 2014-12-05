/** 
 * @file IBackwardGrouper.cpp
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
#include "TransformGroupers/IBackwardGrouper.hpp"

// Project includes
//
#include "TransformConfigurators/TransformStepsMacro.h"

namespace GeoMHDiSCC {

namespace Transform {

   IBackwardGrouper::IBackwardGrouper()
      : split(Splitting::Locations::NONE), 
        mcScalarPacks1D(TransformSteps::Backward<Dimensions::Transform::TRA1D>::SCALAR_VARIABLES), 
        mcGradientPacks1D(TransformSteps::Backward<Dimensions::Transform::TRA1D>::GRAD_VARIABLES), 
        mcVectorPacks1D(TransformSteps::Backward<Dimensions::Transform::TRA1D>::VECTOR_VARIABLES), 
        mcVGradientPacks1D(TransformSteps::Backward<Dimensions::Transform::TRA1D>::VGRAD_VARIABLES), 
        mcCurlPacks1D(TransformSteps::Backward<Dimensions::Transform::TRA1D>::CURL_VARIABLES), 
        mcScalarPacks2D(TransformSteps::Backward<Dimensions::Transform::TRA2D>::SCALAR_VARIABLES),
        mcGradientPacks2D(TransformSteps::Backward<Dimensions::Transform::TRA2D>::GRAD_VARIABLES),
        mcVectorPacks2D(TransformSteps::Backward<Dimensions::Transform::TRA2D>::VECTOR_VARIABLES),
        mcVGradientPacks2D(TransformSteps::Backward<Dimensions::Transform::TRA2D>::VGRAD_VARIABLES),
        mcCurlPacks2D(TransformSteps::Backward<Dimensions::Transform::TRA2D>::CURL_VARIABLES)
   {
   }

   IBackwardGrouper::~IBackwardGrouper()
   {
   }

   ArrayI IBackwardGrouper::namePacks1D(const VariableRequirement& varInfo)
   {
      // Create list of packet sizes
      std::set<int>  list;
      int counter = 0;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for first exchange
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               counter += this->mcScalarPacks1D;
            } else
            {
               counter += this->mcVectorPacks1D;
            }
         }

         // add physical gradient field packs for first exchange
         if(infoIt->second.needPhysicalGradient())
         {
            if(infoIt->second.isScalar())
            {
               counter += this->mcGradientPacks1D;
            } else
            {
               counter += this->mcVGradientPacks1D;
            }
         }

         // add physical curl field packs for first exchange
         if(infoIt->second.needPhysicalCurl())
         {
            if(infoIt->second.isScalar())
            {
               throw Exception("Cannot compute curl of scalar field!");
            } else
            {
               counter += this->mcCurlPacks1D;
            }
         }

         // Add in the combine pack sizes
         if(counter > 0)
         {
            // Add combined pack size to list
            list.insert(counter);
         }

         // Add named pack size
         this->mNamedPacks1D.insert(std::make_pair(infoIt->first, counter));

         // reset counter
         counter = 0;
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
      int counter = 0;

      // loop over all variable information
      VariableRequirement::const_iterator infoIt;
      for(infoIt = varInfo.begin(); infoIt != varInfo.end(); infoIt++)
      {
         // add physical field packs for second exchange
         if(infoIt->second.needPhysical())
         {
            if(infoIt->second.isScalar())
            {
               counter += this->mcScalarPacks2D;
            } else
            {
               counter += this->mcVectorPacks2D;
            }
         }

         // add physical gradient field packs for second exchange
         if(infoIt->second.needPhysicalGradient())
         {
            if(infoIt->second.isScalar())
            {
               counter += this->mcGradientPacks2D;
            } else
            {
               counter += this->mcVGradientPacks2D;
            }
         }

         // add physical curl field packs for second exchange
         if(infoIt->second.needPhysicalCurl())
         {
            if(infoIt->second.isScalar())
            {
               throw Exception("Cannot compute curl of scalar field!");
            } else
            {
               counter += this->mcCurlPacks2D;
            }
         }

         // Add in the combine pack sizes
         if(counter > 0)
         {
            // Add combined pack size to list
            list.insert(counter);
         }

         // Add named pack size
         this->mNamedPacks2D.insert(std::make_pair(infoIt->first, counter));

         // reset counter
         counter = 0;
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
      ArrayI packs = this->namePacks1D(varInfo);

      // Resize packs to single value
      packs.resize(1);
      packs.setConstant(0);

      for(std::map<PhysicalNames::Id, int>::const_iterator it = this->mNamedPacks1D.begin(); it != this->mNamedPacks1D.end(); ++it)
      {
         packs(0) = packs(0) + it->second;
      }

      return packs;
   }

   ArrayI IBackwardGrouper::groupPacks2D(const VariableRequirement& varInfo)
   {  
      // Initialise the number of packs
      ArrayI packs = this->namePacks2D(varInfo);

      // Resize packs to single value
      packs.resize(1);
      packs.setConstant(0);

      for(std::map<PhysicalNames::Id, int>::const_iterator it = this->mNamedPacks2D.begin(); it != this->mNamedPacks2D.end(); ++it)
      {
         packs(0) = packs(0) + it->second;
      }

      return packs;
   }

}
}
