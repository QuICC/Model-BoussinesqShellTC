/** \file SimulationResolution.cpp
 *  \brief Source of the simulation resolution object
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
#include "Base/Resolutions/SimulationResolution.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SimulationResolution::SimulationResolution(const ArrayI& phys, const ArrayI& spec)
   {
      // Initialise 2D dimensions
      this->initDim2D(phys, spec);

      ArrayI extPhys(phys.size() + 1);
      extPhys.segment(0,phys.size()) = phys;
      extPhys(extPhys.size()-1) = this->mDim2D.find(DimensionSpace::PHYSICAL)->second.sum();
      this->mSim.insert(std::make_pair(DimensionSpace::PHYSICAL, extPhys));

      ArrayI extSpec(spec.size() + 1);
      extSpec.segment(0,spec.size()) = spec;
      extSpec(extSpec.size()-1) = this->mDim2D.find(DimensionSpace::SPECTRAL)->second.sum();
      this->mSim.insert(std::make_pair(DimensionSpace::SPECTRAL, extSpec));
   }

   SimulationResolution::~SimulationResolution()
   {
   }

   void SimulationResolution::initDim2D(const ArrayI& phys, const ArrayI& spec)
   {  
      //
      // Create physical space second dimensions
      //
      ArrayI phys2D(phys(2));

      phys2D.setConstant(phys(2));

      this->mDim2D.insert(std::make_pair(DimensionSpace::PHYSICAL, phys2D));

      //
      // Create physical space second dimensions
      //
      ArrayI spec2D(spec(2));

      // If spherical harmonics based scheme
      #if defined EPMPHOENIX_SPATIALSCHEME_WSH || defined EPMPHOENIX_SPATIALSCHEME_TSH || defined EPMPHOENIX_SPATIALSCHEME_FDSH || defined EPMPHOENIX_SPATIALSCHEME_TpSH
         for(int k = 0; k < spec(2); k++)
         {
            spec2D(k) = k + 1;
         }

      // ... else simply set to constant
      #else 
         spec2D.setConstant(spec(2));
      #endif // defined EPMPHOENIX_SPATIALSCHEME_WSH || defined EPMPHOENIX_SPATIALSCHEME_TSH || defined EPMPHOENIX_SPATIALSCHEME_FDSH || defined EPMPHOENIX_SPATIALSCHEME_TpSH

      this->mDim2D.insert(std::make_pair(DimensionSpace::SPECTRAL, spec2D));
   }

   int SimulationResolution::dim(const DimensionSpace::Id id, const int dim) const
   {
      // Safety assertion
      assert(this->mSim.find(id)->second.size() - 1 > dim);

      return this->mSim.find(id)->second(dim);
   }

   int SimulationResolution::dim2D(const DimensionSpace::Id id, const int j) const
   {
      // Safety assertion
      assert(this->mDim2D.find(id)->second.size() > j);

      return this->mDim2D.find(id)->second(j);
   }

   int SimulationResolution::nSlow(const DimensionSpace::Id id) const
   {
      return this->mSim.find(id)->second.tail(1)(0);
   }

   int SimulationResolution::nFast(const DimensionSpace::Id id) const
   {
      return this->mSim.find(id)->second(0);
   }

}
