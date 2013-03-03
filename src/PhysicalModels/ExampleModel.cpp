/** \file ExampleModel.cpp
 *  \brief Source of an example physical model
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

   void ExampleModel::addEquations(SharedSimulation spSim)
   {
      // Add equation
      //pSim->addScalarEquation<AN_EQUATION>();
      
      // Add equation
      //pSim->addVectorEquation<AN_EQUATION>();
   }

   void ExampleModel::setInitialStateFile(SharedSimulation spSim)
   {
      // Set initial state file
      //pSim->setInitialState(A_CONFIGURATION);
   }

   void ExampleModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void ExampleModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Add HDF5 output file
      //pSim->addOutputFile(AN_HDF5FILE);
      
      // Add HDF5 output file
      //pSim->addOutputFile(AN_HDF5FILE);
   }

}
