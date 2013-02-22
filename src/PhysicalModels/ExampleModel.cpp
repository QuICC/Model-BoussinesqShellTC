/** \file ExampleModel.cpp
 *  \brief Source of the dummy simulation for the cylindrical system
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/ExampleModel.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SharedSimulation ExampleModel::createSimulation()
   {
      // Create simulation
      SharedSimulation  spSim;

      return spSim;
   }

   void ExampleModel::equations()
   {
   }

   void ExampleModel::xmlInputs()
   {
   }

   void ExampleModel::hdf5Inputs()
   {
   }

   void ExampleModel::asciiOutputs()
   {
   }

   void ExampleModel::hdf5Outputs()
   {
   }

}
