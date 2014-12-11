/** 
 * @file Resolution.cpp
 * @brief Source of the resolution object for several CPUs
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

   Resolution::Resolution(const std::vector<SharedCoreResolution>& coreRes, const ArrayI& simDim, const ArrayI& transDim)
      : mCores(coreRes)
   {
      // Assert simulation dimensions are the same as cpu dimensions
      assert(simDim.size() == this->cpu()->nDim());

      // Init the simulation resolution
      this->initSimResolution(simDim, transDim);
   }

   Resolution::~Resolution()
   {
   }

   void Resolution::initSimResolution(const ArrayI& simDim, const ArrayI& transDim)
   {
      // Get number of dimensions
      int nDim = this->cpu()->nDim();

      // Create storage for the physical dimensions
      ArrayI phys(nDim);
      for(int i = 0; i < nDim; i++)
      {
         phys(i) = this->cpu()->dim(static_cast<Dimensions::Transform::Id>(i))->dim<Dimensions::Data::DATF1D>();
      }

      ArrayI spec = simDim;
      spec.array() += 1;

      // Create shared pointer of simulation resolution
      this->mspSim = SharedSimulationResolution(new SimulationResolution(phys, spec, transDim));
   }

   int Resolution::nCpu() const
   {
      return this->mCores.size();
   }

   void Resolution::setBoxScale(const Array& boxScale)
   {
      this->mspSim->setBoxScale(boxScale);

      // Set the boxscale on the 
      std::map<Dimensions::Transform::Id,Transform::SharedTransformSetup>::iterator it;
      for(it = this->mTSetups.begin(); it != this->mTSetups.end(); ++it)
      {
         it->second->setBoxScale(this->sim()->boxScale(static_cast<Dimensions::Simulation::Id>(static_cast<int>(it->first))));
      }
      
   }

   void Resolution::setIndexCounter(SharedIndexCounter spCounter)
   {
      this->mspCounter = spCounter;
   }

   SharedIndexCounter   Resolution::counter()
   {
      // Safety assert
      assert(this->mspCounter);

      return this->mspCounter;
   }

   SharedCSimulationResolution Resolution::sim() const
   {
      // Safety assert
      assert(this->mspSim);

      return this->mspSim;
   }

   SharedCCoreResolution Resolution::cpu() const
   {
      // Safety assert
      assert(FrameworkMacro::id() >= 0);
      assert(this->mCores.size() > static_cast<size_t>(FrameworkMacro::id()));

      return this->mCores.at(FrameworkMacro::id());
   }

   SharedCCoreResolution Resolution::cpu(const int id) const
   {
      // Check sizes
      assert(id >= 0);
      assert(static_cast<size_t>(id) < this->mCores.size());
      
      return this->mCores.at(id);
   }

   void Resolution::addTransformSetup(const Dimensions::Transform::Id id, Transform::SharedTransformSetup spSetup)
   {
      this->mTSetups.insert(std::make_pair(id, spSetup));
   }

   Transform::SharedTransformSetup  Resolution::spTransformSetup(const Dimensions::Transform::Id id) const
   {
      // Assert for correct size
      assert(this->mTSetups.size() > static_cast<size_t>(id));

      return this->mTSetups.find(id)->second;
   }

   Datatypes::SharedScalarFieldSetupType Resolution::spFwdSetup(const Dimensions::Transform::Id id) const
   {
      // Iterator for the vector based storages
      std::vector<ArrayI>::const_iterator   vIt;
      
      // Get forward dimensions
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
      spDim1D->setConstant(this->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM));

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

   void Resolution::buildRestriction(std::vector<int>& rSlow, std::vector<std::vector<int> >& rMiddle)
   {
      // Make sure vectors are empty
      rSlow.clear();
      rMiddle.clear();

      #ifdef GEOMHDISCC_MPISPSOLVE
         // Three dimensional matrices
         #if defined GEOMHDISCC_SPATIALSCHEME_TTT
            int k_;
            int j_;
            rSlow.reserve(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
            rMiddle.reserve(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());

            for(int k=0; k < this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               // Extract "physical" index of slow data dimension
               k_ = this->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

               rSlow.push_back(k_);
               rMiddle.push_back(std::vector<int>());
               rMiddle.back().reserve(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k));

               // Loop over middle data dimension
               for(int j=0; j < this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
               {
                  // Extract "physical" index of middle data dimension
                  j_ = this->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);

                  rMiddle.back().push_back(j_);
               }
            }
         // Two dimensional matrices
         #elif defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMHDISCC_SPATIALSCHEME_CFT || defined GEOMHDISCC_SPATIALSCHEME_AFT
            // Make sure middle is empty
            rMiddle.clear();

            rSlow.reserve(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());

            #error "MPI sparse solve restriction for two dimensional matrices not yet implemented."
         #else
            #error "Tried to setup MPI sparse solver for unimplemented scheme."
         #endif //defined GEOMHDISCC_SPATIALSCHEME_TTT
      #else
         throw Exception("buildRestriction was called without MPI sparse solver");
      #endif //GEOMHDISCC_MPISPSOLVE
   }

}
