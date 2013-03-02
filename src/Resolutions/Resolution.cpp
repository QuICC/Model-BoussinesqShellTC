/** \file Resolution.cpp
 *  \brief Source of the resolution object for several CPUs
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <cassert>

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
         phys(i) = this->cpu()->dim(static_cast<Dimensions::Transform::Id>((nDim-1)-i))->dim<Dimensions::Data::DATF1D>();
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

   void Resolution::addTransformSetup(const Dimensions::Transform::Id id, Transform::SharedFftSetup spSetup)
   {
      this->mTSetups.insert(std::make_pair(id, spSetup));
   }

   Transform::SharedFftSetup  Resolution::spTransformSetup(const Dimensions::Transform::Id id) const
   {
      // Assert for correct size
      assert(this->mTSetups.size() > static_cast<size_t>(id));

      return this->mTSetups.find(id)->second;
   }

   Datatypes::SharedScalarFieldSetupType Resolution::spFwdSetup(const Dimensions::Transform::Id id) const
   {
      // Iterator for the vector based storages
      std::vector<ArrayI>::const_iterator   vIt;
      
      // Get backward dimensions
      SharedArrayI   spDim1D(new ArrayI(this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>()));
      for(int i = 0; i < spDim1D->size(); ++i)
      {
         (*spDim1D)(i) = this->cpu()->dim(id)->dim<Dimensions::Data::DATF1D>(i);
      }

      // Get 2D dimensions
      SharedArrayI   spDim2D(new ArrayI(this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>()));
      for(int i = 0; i < spDim2D->size(); ++i)
      {
         (*spDim2D)(i) = this->cpu()->dim(id)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Datatypes::SharedScalarFieldSetupType(new Datatypes::ScalarFieldSetupType(spDim1D, spDim2D, this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>()));
   }

   Datatypes::SharedScalarFieldSetupType Resolution::spBwdSetup(const Dimensions::Transform::Id id) const
   {
      // Iterator for the vector based storages
      std::vector<ArrayI>::const_iterator   vIt;
      
      // Get backward dimensions
      SharedArrayI   spDim1D(new ArrayI(this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>()));
      for(int i = 0; i < spDim1D->size(); ++i)
      {
         (*spDim1D)(i) = this->cpu()->dim(id)->dim<Dimensions::Data::DATB1D>(i);
      }

      // Get 2D dimensions
      SharedArrayI   spDim2D(new ArrayI(this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>()));
      for(int i = 0; i < spDim2D->size(); ++i)
      {
         (*spDim2D)(i) = this->cpu()->dim(id)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Datatypes::SharedScalarFieldSetupType(new Datatypes::ScalarFieldSetupType(spDim1D, spDim2D, this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>()));
   }

   Datatypes::SharedScalarFieldSetupType Resolution::spSpectralSetup() const
   {
      // Iterator for the vector based storages
      std::vector<ArrayI>::const_iterator   vIt;
      
      // Get backward dimensions
      SharedArrayI   spDim1D(new ArrayI(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>()));
      spDim1D->setConstant(this->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));

      // Get 2D dimensions
      SharedArrayI   spDim2D(new ArrayI(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>()));
      for(int i = 0; i < spDim2D->size(); ++i)
      {
         (*spDim2D)(i) = this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return Datatypes::SharedScalarFieldSetupType(new Datatypes::ScalarFieldSetupType(spDim1D, spDim2D, this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>()));
   }

   Datatypes::SharedScalarFieldSetupType Resolution::spPhysicalSetup() const
   {
      return this->spFwdSetup(static_cast<Dimensions::Transform::Id>(this->cpu()->nDim()-1));
   }

}
