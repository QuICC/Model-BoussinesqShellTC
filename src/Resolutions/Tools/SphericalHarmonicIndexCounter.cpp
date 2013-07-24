/** \file SphericalHarmonicIndexCounter.cpp
 *  \brief Source of spherical harmonic index counter
 */

// System includes
//

// External includes
//

// Class include
//
#include "Resolutions/Tools/SphericalHarmonicIndexCounter.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SphericalHarmonicIndexCounter::SphericalHarmonicIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu)
      : IndexCounter(), mspSim(spSim), mspCpu(spCpu)
   {
   }

   SphericalHarmonicIndexCounter::~SphericalHarmonicIndexCounter()
   {
   }

   ArrayI SphericalHarmonicIndexCounter::orderedDimensions(const Dimensions::Space::Id spaceId) const
   {
      ArrayI dims = this->mspSim->dimensions(spaceId);

      return this->orderedDimensions(dims, spaceId);
   }

   ArrayI SphericalHarmonicIndexCounter::orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const
   {
      // Storage for the ordered dimensions
      ArrayI oDims(dims.size());

      // In spectral space reduce dimensions to 1D, NH (=number of harmonics)
      if(spaceId == Dimensions::Space::SPECTRAL)
      {
         assert(dims.size() == 3);

         oDims.resize(2);
         oDims(0) = dims(0);
         oDims(1) = 0;

         for(int l = 0; l < dims(1); ++l)
         {
            for(int m = 0; m < std::min(l+1,dims(2)); ++m)
            {
               oDims(1)++;
            }
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

   void SphericalHarmonicIndexCounter::computeOffsets(std::vector<std::vector<SphericalHarmonicIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const
   {
      this->computeOffsets(offsets, spaceId, this->mspSim);
   }

   void SphericalHarmonicIndexCounter::computeOffsets(std::vector<std::vector<SphericalHarmonicIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const
   {
      Dimensions::Transform::Id transId;
      Dimensions::Simulation::Id simId;
      
      // Clear the vector of offsets
      offsets.clear();
      std::vector<OffsetType>  offV;

      // In spectral space offset computation, spherical harmonic triangular truncation make it complicated
      if(spaceId == Dimensions::Space::SPECTRAL)
      {
         transId = Dimensions::Transform::TRA1D;
         simId = Dimensions::Simulation::SIM2D;

         // Loop over all local harmonic degrees l 
         OffsetType offset = 0;
         int l0 = 0;
         for(int iL = 0; iL < this->mspCpu->dim(transId)->dim<Dimensions::Data::DAT3D>(); ++iL)
         {
            int l_ = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT3D>(iL);
            if(l_ < spRef->dim(simId,spaceId))
            {
               // Compute the offset to the local harmonic degree l - 1
               for(int l = l0; l < l_; ++l)
               {  
                  for(int m = 0; m < std::min(l+1,spRef->dim(Dimensions::Simulation::SIM3D,spaceId)); ++m)
                  {
                     offset++;
                  }  
               }

               // Compute offset for the local m
               offV.clear();
               for(int iM = 0; iM < this->mspCpu->dim(transId)->dim<Dimensions::Data::DAT2D>(iL); ++iM)
               {
                  int m_ = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT2D>(iM, iL);
                  if(m_ < spRef->dim(Dimensions::Simulation::SIM3D,spaceId))
                  {
                     offV.push_back(offset + m_);
                  }
               }

               offsets.push_back(offV);

               l0 = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT3D>(iL);
            }
         }

      //  Physical space offset computation (regular)
      } else //if(spaceId == Dimensions::Space::PHYSICAL)
      {
         transId = Dimensions::Transform::TRA3D;
         simId = Dimensions::Simulation::SIM1D;

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

               offsets.push_back(offV);
            }
         }
      }
   }
}