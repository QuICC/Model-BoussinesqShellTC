/**
 * @file StateGenerator.hpp
 * @brief High level implementation of a general state generator 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STATEGENERATOR_HPP
#define STATEGENERATOR_HPP

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
    * @brief High level implementation of a general state generator
    */
   class StateGenerator: public SimulationBase
   {
      public:
         /**
          * @brief Constructor
          */
         StateGenerator();

         /**
          * @brief Simple empty destructor
          */
         virtual ~StateGenerator();

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

   /// Typedef for a shared pointer of a StateGenerator
   typedef SharedPtrMacro<StateGenerator> SharedStateGenerator;

}

#endif // STATEGENERATOR_HPP
