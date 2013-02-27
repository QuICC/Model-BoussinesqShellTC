/** \file Resolution.cpp
 *  \brief Source of the resolution object for several CPUs
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "Resolutions/Resolution.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   Resolution::Resolution(const std::vector<SharedCoreResolution>& coreRes, const ArrayI& simDim)
      : mCores(coreRes)
   {
      // Assert simulation dimensions are the same as cpu dimensions
      assert(simDim.size() == this->cpu()->nDim());

      // Init the simulation resolution
      this->initSimResolution(simDim);
   }

   Resolution::~Resolution()
   {
   }

   void Resolution::initSimResolution(const ArrayI& simDim)
   {
      // Get number of dimensions
      int nDim = this->cpu()->nDim();

      // Create storage for the physical dimensions
      ArrayI phys(nDim);
      for(int i = 0; i < nDim; i++)
      {
         phys(i) = this->cpu()->dim(static_cast<Dimensions::Transform::Id>((nDim-1)-i))->dimFwd();
      }

      ArrayI spec = simDim;
      spec.array() += 1;

      // Create shared pointer of simulation resolution
      this->mspSim = SharedSimulationResolution(new SimulationResolution(phys, spec));
   }

   int Resolution::nCpu() const
   {
      return this->mCores.size();
   }

   SharedCSimulationResolution Resolution::sim() const
   {
      return this->mspSim;
   }

   SharedCCoreResolution Resolution::cpu() const
   {
      // Safety assert
      assert(FrameworkMacro::id() >= 0);
      assert(this->mCores.size() > static_cast<unsigned int>(FrameworkMacro::id()));

      return this->mCores.at(FrameworkMacro::id());
   }

   SharedCCoreResolution Resolution::cpu(const int id) const
   {
      // Check sizes
      assert(id >= 0);
      assert(static_cast<unsigned int>(id) < this->mCores.size());
      
      return this->mCores.at(id);
   }

   Datatypes::SharedScalarFieldSetupType Resolution::spFwdSetup() const
   {
      return this->cpu()->dim(static_cast<Dimensions::Transform::Id>(this->cpu()->nDim()-1))->spFwdSetup();
   }

   Datatypes::SharedScalarFieldSetupType Resolution::spFwdSetup(const Dimensions::Transform::Id id) const
   {
      return this->cpu()->dim(id)->spFwdSetup();
   }

   Datatypes::SharedScalarFieldSetupType Resolution::spBwdSetup() const
   {
      return this->cpu()->dim(Dimensions::Transform::TRA1D)->spBwdSetup();
   }

   Datatypes::SharedScalarFieldSetupType Resolution::spBwdSetup(const Dimensions::Transform::Id id) const
   {
      return this->cpu()->dim(id)->spBwdSetup();
   }

}
