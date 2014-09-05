/** 
 * @file DoublePeriodicIndexCounter.cpp
 * @brief Source of double periodic index counter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Resolutions/Tools/DoublePeriodicIndexCounter.hpp"

// Project includes
//

#include <iostream>
namespace GeoMHDiSCC {

   DoublePeriodicIndexCounter::DoublePeriodicIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu)
      : IndexCounter(), mspSim(spSim), mspCpu(spCpu)
   {
   }

   DoublePeriodicIndexCounter::~DoublePeriodicIndexCounter()
   {
   }

   ArrayI DoublePeriodicIndexCounter::orderedDimensions(const Dimensions::Space::Id spaceId) const
   {
      ArrayI dims = this->mspSim->dimensions(spaceId);

      return this->orderedDimensions(dims, spaceId);
   }

   ArrayI DoublePeriodicIndexCounter::orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const
   {
      // Storage for the ordered dimensions
      ArrayI oDims(dims.size());

      // Spectral and transform space ordering is 1D, 3D, 2D
      if(spaceId == Dimensions::Space::SPECTRAL || spaceId == Dimensions::Space::TRANSFORM)
      {
         oDims(0) = dims(0);
         for(int i = 1; i < dims.size(); ++i)
         {
            oDims(i) = dims(dims.size()-i);
         }
      
      //  Physical space ordering is 3D, 2D, 1D
      } else //if(spaceId == Dimensions::Space::PHYSICAL)
      {
         for(int i = 0; i < dims.size(); ++i)
         {
            oDims(i) = dims(dims.size()-1-i);
         }
      }

      return oDims;
   }

   void DoublePeriodicIndexCounter::computeOffsets(std::vector<std::vector<DoublePeriodicIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const
   {
      this->computeOffsets(offsets, spaceId, this->mspSim);
   }

   void DoublePeriodicIndexCounter::computeOffsets(std::vector<std::vector<DoublePeriodicIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const
   {
      Dimensions::Transform::Id transId;
      Dimensions::Simulation::Id simId;

      // Select transform dimension depending on dimension space
      if(spaceId == Dimensions::Space::SPECTRAL || spaceId == Dimensions::Space::TRANSFORM)
      {
         transId = Dimensions::Transform::TRA1D;
         simId = Dimensions::Simulation::SIM2D;

      } else //if(spaceId == Dimensions::Space::PHYSICAL)
      {
         transId = Dimensions::Transform::TRA3D;
         simId = Dimensions::Simulation::SIM1D;
      }

      // Clear the vector of offsets
      offsets.clear();
      std::vector<OffsetType>  offV;

      offV.push_back(0);
      offV.push_back(0);

      // Select transform dimension depending on dimension space
      if(spaceId == Dimensions::Space::SPECTRAL)
      {
         // Get full slowest resolution resolution
         int dat3D = this->mspSim->dim(simId, Dimensions::Space::TRANSFORM);
         int ref3D = spRef->dim(simId,spaceId)/2 + 1;

         for(int i=0; i < this->mspCpu->dim(transId)->dim<Dimensions::Data::DAT3D>(); ++i)
         {
            int i_ = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT3D>(i);
            // Check if value is available in file
            if(i_ < ref3D)
            {
               // Compute offset for third dimension
               offV.at(0) = i_;

               // Compute offset for second dimension
               offV.at(1) = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT2D>(0,i);

               offsets.push_back(offV);
               std::cerr << i_ << std::endl;
            } else if(dat3D - i_ < ref3D)
            {
               // Compute offset for third dimension
               offV.at(0) = i_ - (dat3D - spRef->dim(simId,spaceId));

               // Compute offset for second dimension
               offV.at(1) = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT2D>(0,i);

               offsets.push_back(offV);
               std::cerr << i_ - (dat3D - spRef->dim(simId,spaceId)) << std::endl;
            }
         }

      } else //if(spaceId == Dimensions::Space::PHYSICAL || spaceId == Dimensions::Space::TRANSFORM)
      {
         for(int i=0; i < this->mspCpu->dim(transId)->dim<Dimensions::Data::DAT3D>(); ++i)
         {
            int i_ = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT3D>(i);
            // Check if value is available in file
            if(i_ < spRef->dim(simId,spaceId))
            {
               // Compute offset for third dimension
               offV.at(0) = i_;

               // Compute offset for second dimension
               offV.at(1) = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT2D>(0,i);

               offsets.push_back(offV);
            }
         }
      }
   }
}
