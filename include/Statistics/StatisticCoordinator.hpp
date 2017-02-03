/**
 * @file StatisticCoordinator.hpp
 * @brief Coordinator for the statistics computations 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STATISTICCOORDINATOR_HPP
#define STATISTICCOORDINATOR_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <vector>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace QuICC {

namespace Statistics {

   /**
    * @brief Coordinator for the statistics computations
    */
   class StatisticCoordinator
   {
      public:
         /**
          * @brief Constructor
          */
         StatisticCoordinator();

         /**
          * @brief Constructor
          */
         ~StatisticCoordinator();

         /**
          * @brief Initialise the coordinator
          */
         void init(); 

      protected:

      private:
   };
}
}

#endif // STATISTICCOORDINATOR_HPP
