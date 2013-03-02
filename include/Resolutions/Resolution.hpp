/** \file Resolution.hpp
 *  \brief Definition of a resolution object
 *
 *  \mhdBug Needs test
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
#include "FastTransforms/FftSetup.hpp"

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
          * @param coreRes Resolution object for the different CPUs
          * @param simDim  Simulation dimensions
          */
         Resolution(const std::vector<SharedCoreResolution>& coreRes, const ArrayI& simDim);

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
         Transform::SharedFftSetup spTransformSetup(const Dimensions::Transform::Id dim) const;

         /**
          * @brief Add a transform setup 
          *
          * @param dim     Dimension corresponding to setup
          * @param spSetup Transform setup
          */
         void addTransformSetup(const Dimensions::Transform::Id dim, Transform::SharedFftSetup spSetup);

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

      protected:

      private:
         /**
          * @brief Initialise the simulation resolution
          *
          * @param simDim  Simulation dimensions
          */
         void initSimResolution(const ArrayI& simDim);

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
         std::map<Dimensions::Transform::Id,Transform::SharedFftSetup>   mTSetups;
   };

   /// Typedef for a shared pointer to a Resolution object
   typedef SharedPtrMacro<Resolution>   SharedResolution;

   /// Typedef for a shared pointer to a const Resolution object
   typedef SharedPtrMacro<const Resolution>   SharedCResolution;

}

#endif // RESOLUTION_HPP
