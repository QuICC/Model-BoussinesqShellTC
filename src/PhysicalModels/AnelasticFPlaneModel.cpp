/** 
 * @file AnelasticFPlaneModel.cpp
 * @brief Source of the anelastic F-plane physical model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalModels/AnelasticFPlaneModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/FPlane/Anelastic/AnelasticFPlaneStreamfunction.hpp"
#include "Equations/FPlane/Anelastic/AnelasticFPlaneVertical.hpp"
#include "Equations/FPlane/Anelastic/AnelasticFPlaneTransport.hpp"
#include "PhysicalModels/PhysicalModelBase.hpp"

namespace GeoMHDiSCC {

   const std::string AnelasticFPlaneModel::PYMODULE = "anelastic_fplane_model";

   const std::string AnelasticFPlaneModel::PYCLASS = "AnelasticFPlaneModel";

   void AnelasticFPlaneModel::addEquations(SharedSimulation spSim)
   {
      throw Exception("Not implemented yet!");
   }

   void AnelasticFPlaneModel::addStates(SharedStateGenerator spGen)
   {
      throw Exception("Not implemented yet!");
   }

   void AnelasticFPlaneModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      throw Exception("Not implemented yet!");
   }

   void AnelasticFPlaneModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      throw Exception("Not implemented yet!");
   }

   void AnelasticFPlaneModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void AnelasticFPlaneModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addHdf5OutputFile(spState);
   }

   void AnelasticFPlaneModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = PhysicalModelBase::fieldIds();

      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spInit(new IoVariable::StateFileReader("_initial", SchemeType::type(), SchemeType::isRegular()));

      // Set expected field names
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spInit->expect(*it);
      }

      // Set simulation state
      spSim->setInitialState(spInit);
   }

}
