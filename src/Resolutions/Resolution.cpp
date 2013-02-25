/** \file Resolution.cpp
 *  \brief Source of the resolution object for several CPUs
 */

// Configuration includes
//
#include "PrepMacros/FrameworkMacro.h"

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

   Resolution::Resolution(const std::vector<SharedCCoreResolution>& coreRes, const ArrayI& simDim)
      : mCores(coreRes)
   {
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
         phys(i) = this->cpu()->dim((nDim-1)-i)->dimFwd();
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
      assert(this->mCores.size() > FrameworkMacro::id());

      return this->mCores.at(FrameworkMacro::id());
   }

   SharedCCoreResolution Resolution::cpu(const int id) const
   {
      // Check sizes
      assert(id < this->mCores.size());
      
      return this->mCores.at(id);
   }

   Datatypes::SharedScalarFieldSetup Resolution::spFwdSetup() const
   {
      return this->cpu()->dim(this->cpu()->nDim()-1)->spFwdSetup();
   }

   Datatypes::SharedScalarFieldSetup Resolution::spFwdSetup(const int dim) const
   {
      return this->cpu()->dim(dim)->spFwdSetup();
   }

   Datatypes::SharedScalarFieldSetup Resolution::spBwdSetup() const
   {
      return this->cpu()->dim(0)->spBwdSetup();
   }

   Datatypes::SharedScalarFieldSetup Resolution::spBwdSetup(const int dim) const
   {
      return this->cpu()->dim(dim)->spBwdSetup();
   }

}
