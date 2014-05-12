/** 
 * @file BoussinesqAnnulusModel.cpp
 * @brief Source of the Boussinesq annulus physical model
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
#include "PhysicalModels/BoussinesqAnnulusModel.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "IoVariable/StateFileReader.hpp"
#include "IoVariable/StateFileWriter.hpp"
#include "IoVariable/VisualizationFileWriter.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqAnnulusTransport.hpp"
#include "Equations/Annulus/Boussinesq/BoussinesqAnnulusVelocity.hpp"
#include "Generator/States/AnnulusExactScalarState.hpp"
#include "Generator/States/AnnulusExactVectorState.hpp"
#include "Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "Generator/Visualizers/VectorFieldVisualizer.hpp"

namespace GeoMHDiSCC {

   const std::string BoussinesqAnnulusModel::PYNAME = "boussinesq_annulus_model";

   std::vector<PhysicalNames::Id> BoussinesqAnnulusModel::fieldIds()
   {
      // Create storage
      std::vector<PhysicalNames::Id> ids;

      // Add temperature
      ids.push_back(PhysicalNames::TEMPERATURE);

      // Add temperature
      ids.push_back(PhysicalNames::VELOCITY);

      return ids;
   }

   std::vector<NonDimensional::Id> BoussinesqAnnulusModel::paramIds()
   {
      // Create storage
      std::vector<NonDimensional::Id> ids;

      // Add gap width
      ids.push_back(NonDimensional::GAPWIDTH);

      // Add radii ratio
      ids.push_back(NonDimensional::RRATIO);

      return ids;
   }

   std::vector<bool> BoussinesqAnnulusModel::isPeriodicBox()
   {
      std::vector<bool> box;

      // X direction is not periodic box
      box.push_back(false);

      // Y direction is periodic box
      box.push_back(false);

      // Z direction is not periodic box
      box.push_back(false);

      return box;
   }

   void BoussinesqAnnulusModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::BoussinesqAnnulusTransport>(BoussinesqAnnulusModel::PYNAME);
      
      // Add Navier-Stokes equation
      spSim->addVectorEquation<Equations::BoussinesqAnnulusVelocity>(BoussinesqAnnulusModel::PYNAME);
   }

   void BoussinesqAnnulusModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedAnnulusExactScalarState spSExact;
      Equations::SharedAnnulusExactVectorState spVExact;

      // Add temperature initial state generator
      spSExact = spGen->addScalarEquation<Equations::AnnulusExactScalarState>();
      spSExact->setIdentity(PhysicalNames::TEMPERATURE);
      spSExact->setStateType(Equations::AnnulusExactScalarState::CONSTANT);

      // Add temperature initial state generator
      spVExact = spGen->addVectorEquation<Equations::AnnulusExactVectorState>();
      spVExact->setIdentity(PhysicalNames::VELOCITY);
      spVExact->setStateType(Equations::AnnulusExactVectorState::CONSTANT);

      // Add output file
      IoVariable::SharedStateFileWriter spOut(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spOut->expect(PhysicalNames::VELOCITY);
      spGen->addOutputFile(spOut);
   }

   void BoussinesqAnnulusModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spField;

      // Add first field visualization
      spField = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spField->setFields(true, false);
      spField->setIdentity(PhysicalNames::TEMPERATURE);

      // Add output file
      IoVariable::SharedVisualizationFileWriter spOut(new IoVariable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::TEMPERATURE);
      spVis->addOutputFile(spOut);
   }

   void BoussinesqAnnulusModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      IoVariable::SharedStateFileReader spIn(new IoVariable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::TEMPERATURE);

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void BoussinesqAnnulusModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
      
      // Add ASCII output file
      //pSim->addOutputFile(AN_ASCIIFILE);
   }

   void BoussinesqAnnulusModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqAnnulusModel::fieldIds();

      // Create and add state file to IO
      IoVariable::SharedStateFileWriter spState(new IoVariable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }
      spSim->addOutputFile(spState);
   }

   void BoussinesqAnnulusModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<GeoMHDiSCC::PhysicalNames::Id>::const_iterator  it;
      std::vector<GeoMHDiSCC::PhysicalNames::Id> ids = BoussinesqAnnulusModel::fieldIds();

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
