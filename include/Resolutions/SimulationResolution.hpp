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
         virtual ~SimulationResolution();

         /**
          * @brief Get simulation's dimensions
          *
          * @param id ID of the space (Physical, Spectral)
          * @param dim Dimension id
          */
         int dim(const Dimensions::Space::Id id, const int dim) const;

         /**
          * @brief Get simulation's second dimension
          *
          * @param id   ID of the space (Physical, Spectral)
          * @param j    Index of the second dimension
          */
         int dim2D(const Dimensions::Space::Id id, const int j) const;

      protected:

      private:
         /**
          * brief Storage
          */
         std::map<Dimensions::Space::Id, ArrayI> mSim;

         /**
          * brief Storage
          */
         std::map<Dimensions::Space::Id, ArrayI> mDim2D;

         /**
          * @brief Create second dimension values
          *
          * @param phys    Dimensions for the physical space
          * @param spec    Dimensions for the spectral space
          */
         void initDim2D(const ArrayI& phys, const ArrayI& spec);
   };


   /// Typedef for a shared pointer to a SimulationResolution object
   typedef SharedPtrMacro<SimulationResolution>   SharedSimulationResolution;

   /// Typedef for a shared pointer to a const SimulationResolution object
   typedef SharedPtrMacro<const SimulationResolution>   SharedCSimulationResolution;

}

#endif // SIMULATIONRESOLUTION_HPP
