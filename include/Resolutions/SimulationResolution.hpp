/** \file SimulationResolution.hpp
 *  \brief Definition of a simulation resolution object
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

         /**
          * @brief Get the simulation dimensions ordered by index access speed (fast -> slow)
          */
         ArrayI orderedDims(const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Get the box scale
          *
          * @para id ID of the simulation dimension
          */
         MHDFloat boxScale(const Dimensions::Simulation::Id id) const;

         /**
          * @brief Set the box scale (if not it will be initialised to 1)
          */
         void setBoxScale(const Array& boxScale);

      protected:

      private:
         /**
          * brief Storage
          */
         std::map<Dimensions::Space::Id, ArrayI> mDim;

         /**
          * @brief Storage for the box scale
          */
         Array mBoxScale;
   };


   /// Typedef for a shared pointer to a SimulationResolution object
   typedef SharedPtrMacro<SimulationResolution>   SharedSimulationResolution;

   /// Typedef for a shared pointer to a const SimulationResolution object
   typedef SharedPtrMacro<const SimulationResolution>   SharedCSimulationResolution;

}

#endif // SIMULATIONRESOLUTION_HPP
