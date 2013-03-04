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
#include "Enums/PhysicalNames.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"

namespace GeoMHDiSCC {

   void ExampleModel::addEquations(SharedSimulation spSim)
   {
      // Add equation
      //pSim->addScalarEquation<AN_EQUATION>();
      
      // Add equation
      //pSim->addVectorEquation<AN_EQUATION>();
   }

   void ExampleModel::setInitialState(SharedSimulation spSim)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spInit(new IoVariable::StateFileReader("_initial", SchemeType::type()));
      spInit->expect(PhysicalNames::STREAMFUNCTION);
      spInit->expect(PhysicalNames::VELOCITYZ);
      spInit->expect(PhysicalNames::TEMPERATURE);
      spSim->setInitialState(spInit);
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
      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type()));
      spState->expect(PhysicalNames::STREAMFUNCTION);
      spState->expect(PhysicalNames::VELOCITYZ);
      spState->expect(PhysicalNames::TEMPERATURE);
      spSim->addOutputFile(spState);
      
      // Create and add visualization file to IO
      IoVariable::SharedVisualizationFileWriter spViz(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spViz->expect(PhysicalNames::STREAMFUNCTION);
      spViz->expect(PhysicalNames::VELOCITYZ);
      spViz->expect(PhysicalNames::TEMPERATURE);
      spSim->addOutputFile(spViz);
   }

}
