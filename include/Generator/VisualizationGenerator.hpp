/** \file VisualizationGenerator.hpp
 *  \brief High level implementation of a general visualization file generator
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
          * @brief Allow for implementation specific output tuning
          */
         virtual void tuneOutput();

         /**
          * @brief Write the output if required
          */
         void writeOutput();
   };

   /// Typedef for a shared pointer of a VisualizationGenerator
   typedef SharedPtrMacro<VisualizationGenerator> SharedVisualizationGenerator;

}

#endif // VISUALIZATIONGENERATOR_HPP
