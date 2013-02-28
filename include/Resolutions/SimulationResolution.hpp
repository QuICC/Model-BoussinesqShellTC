/** \file SimulationResolution.hpp
 *  \brief Definition of a simulation resolution object
 *
 *  \mhdBug Needs test
 */

#ifndef SIMULATIONRESOLUTION_HPP
#define SIMULATIONRESOLUTION_HPP

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
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Definition of a simulation resolution object
    */
   class SimulationResolution
   {
      public:
         /**
          * @brief Constructor
          *
          * @param phys    Dimensions for the physical space
          * @param spec    Dimensions for the spectral space
          */
         SimulationResolution(const ArrayI& phys, const ArrayI& spec);

         /**
          * @brief Empty Destructor
          */
         ~SimulationResolution();

         /**
          * @brief Get simulation's dimensions
          *
          * @param simId  ID of the simulation dimension (SIM1D, SIM2D, SIM3D)
          * @param spaceId ID of the space (PHYSICAL, SPECTRAL)
          */
         int dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId) const;

      protected:

      private:
         /**
          * brief Storage
          */
         std::map<Dimensions::Space::Id, ArrayI> mDim;
   };


   /// Typedef for a shared pointer to a SimulationResolution object
   typedef SharedPtrMacro<SimulationResolution>   SharedSimulationResolution;

   /// Typedef for a shared pointer to a const SimulationResolution object
   typedef SharedPtrMacro<const SimulationResolution>   SharedCSimulationResolution;

}

#endif // SIMULATIONRESOLUTION_HPP
