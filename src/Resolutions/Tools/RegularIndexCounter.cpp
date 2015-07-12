/** 
 * @file RegularIndexCounter.cpp
 * @brief Source of regular index counter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Resolutions/Tools/RegularIndexCounter.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   RegularIndexCounter::RegularIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu)
      : IndexCounter(), mspSim(spSim), mspCpu(spCpu)
   {
   }

   RegularIndexCounter::~RegularIndexCounter()
   {
   }

   ArrayI RegularIndexCounter::orderedDimensions(const Dimensions::Space::Id spaceId) const
   {
      ArrayI dims = this->mspSim->dimensions(spaceId);

      return this->orderedDimensions(dims, spaceId);
   }

   ArrayI RegularIndexCounter::orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const
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

   void RegularIndexCounter::computeOffsets(std::vector<std::vector<RegularIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const
   {
      this->computeOffsets(offsets, spaceId, this->mspSim);
   }

   void RegularIndexCounter::computeOffsets(std::vector<std::vector<RegularIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const
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
         transId = Dimensions::Transform::TRAND;
         simId = Dimensions::Simulation::SIM1D;
      }

      // Clear the vector of offsets
      offsets.clear();
      std::vector<OffsetType>  offV;

      offV.push_back(0);
      offV.push_back(0);
      offV.push_back(0);
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

            // Store 3D index
            offV.at(2) = i;

            offsets.push_back(offV);
         }
      }
   }
}
