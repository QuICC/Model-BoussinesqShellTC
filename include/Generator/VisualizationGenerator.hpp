/**
 * @file VisualizationGenerator.hpp
 * @brief High level implementation of a general visualization file generator 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VISUALIZATIONGENERATOR_HPP
#define VISUALIZATIONGENERATOR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Simulation/SimulationBase.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief High level implementation of a general visualization file generator
    */
   class VisualizationGenerator: public SimulationBase
   {
      public:
         /**
          * @brief Constructor
          */
         VisualizationGenerator();

         /**
          * @brief Simple empty destructor
          */
         virtual ~VisualizationGenerator();

      protected:

      private:
         /**
          * @brief Do operations required just before starting the main loop
          */
         virtual void preRun();

         /**
          * @brief Do operations required during the main loop
          */
         virtual void mainRun();

         /**
          * @brief Do operations required just after finishing the main loop
          */
         virtual void postRun();

         /**
          * @brief Write the output if required
          */
         void writeOutput();

         /**
          * @brief Allow for implementation specific output tuning
          */
         virtual void tuneOutput();

         /**
          * @brief Allow for additional operators on the initial state input file
          *
          * @param spInitFile Shared initial state file
          */
         virtual void tuneInitialState(IoVariable::SharedStateFileReader spInitFile);
   };

   /// Typedef for a shared pointer of a VisualizationGenerator
   typedef SharedPtrMacro<VisualizationGenerator> SharedVisualizationGenerator;

}

#endif // VISUALIZATIONGENERATOR_HPP
