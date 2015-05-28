/** 
 * @file Resolution.hpp
 * @brief Definition of a resolution object
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef RESOLUTION_HPP
#define RESOLUTION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Resolutions/SimulationResolution.hpp"
#include "Resolutions/CoreResolution.hpp"
#include "Resolutions/TransformSetup.hpp"
#include "Resolutions/Tools/IndexCounter.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Definition of a resolution object
    */
   class Resolution
   {
      public:
         /**
          * @brief Constructor
          *
          * @param coreRes    Resolution object for the different CPUs
          * @param simDim     Simulation dimensions
          * @param transDim   Transform dimensions
          */
         Resolution(const std::vector<SharedCoreResolution>& coreRes, const ArrayI& simDim, const ArrayI& transDim);

         /**
          * @brief Empty Destructor
          */
         ~Resolution();

         /**
          * @brief Get the simulation resolution
          */
         SharedCSimulationResolution sim() const;

         /**
          * @brief Get resolution corresponding to local CPU
          */
         SharedCCoreResolution cpu() const;

         /**
          * @brief Get resolution corresponding to id
          *
          * @param id CPU id
          */
         SharedCCoreResolution cpu(const int id) const;

         /**
          * @brief Get Number of CPUs
          */
         int nCpu() const;

         /**
          * @brief Get the transform setup 
          *
          * @param dim Dimension for which to get setup
          */
         Transform::SharedTransformSetup spTransformSetup(const Dimensions::Transform::Id dim) const;

         /**
          * @brief Add a transform setup 
          *
          * @param dim     Dimension corresponding to setup
          * @param spSetup Transform setup
          */
         void addTransformSetup(const Dimensions::Transform::Id dim, Transform::SharedTransformSetup spSetup);

         /**
          * @brief Get the forward scalar field setup
          *
          * @param dim Dimension for which to get setup
          */
         Datatypes::SharedScalarFieldSetupType spFwdSetup(const Dimensions::Transform::Id dim) const;

         /**
          * @brief Get the backward scalar field setup
          *
          * @param dim Dimension for which to get setup
          */
         Datatypes::SharedScalarFieldSetupType spBwdSetup(const Dimensions::Transform::Id dim) const;

         /**
          * @brief Get the forward scalar field setup for last dimension
          */
         Datatypes::SharedScalarFieldSetupType spPhysicalSetup() const;

         /**
          * @brief Get the backward scalar field setup for the first dimension
          */
         Datatypes::SharedScalarFieldSetupType spSpectralSetup() const;

         /**
          * @brief Set the box scale for the periodic box dimensions
          */
         void setBoxScale(const Array& boxScale);

         /**
          * @brief Set the index counter
          */
         void setIndexCounter(SharedIndexCounter spCounter);

         /**
          * @brief Get the index counter
          */
         SharedIndexCounter counter();

         /**
          * @brief Build restriction for MPI sparse solver
          */
         void buildRestriction(std::vector<int>& rSlow, std::vector<std::vector<int> >& rMiddle, const int k);

      protected:

      private:
         /**
          * @brief Initialise the simulation resolution
          *
          * @param simDim     Simulation dimensions
          * @param transDim   Transform dimensions
          */
         void initSimResolution(const ArrayI& simDim, const ArrayI& transDim);

         /**
          * @brief Shared simulation resolution
          */
         SharedSimulationResolution mspSim;

         /**
          * @brief Storage for all the core resolutions
          */
         std::vector<SharedCoreResolution>   mCores;

         /**
          * @brief Storage for the transform setups
          */
         std::map<Dimensions::Transform::Id,Transform::SharedTransformSetup>   mTSetups;

         /**
          * @brief Storage for the index counter
          */
         SharedIndexCounter   mspCounter;
   };

   /// Typedef for a shared pointer to a Resolution object
   typedef SharedPtrMacro<Resolution>   SharedResolution;

   /// Typedef for a shared pointer to a const Resolution object
   typedef SharedPtrMacro<const Resolution>   SharedCResolution;

}

#endif // RESOLUTION_HPP
