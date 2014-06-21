/** 
 * @file PhysicalModelBase.hpp
 * @brief Implementation of a base for all physical models
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PHYSICALMODELBASE_HPP
#define PHYSICALMODELBASE_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Simulation/Simulation.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a base for all physical models
    */
   class PhysicalModelBase
   {
      public:
         /**
          * @brief Get vector of names for the boundary conditions
          */
         static std::vector<PhysicalNames::Id> fieldIds();

         /**
          * @brief Get vector of names for the nondimensional parameters
          */
         static std::vector<NonDimensional::Id> paramIds();

         /**
          * @brief Get vector of bools about periodic box
          */
         static std::vector<bool> isPeriodicBox();

      protected:

      private:
         /**
          * @brief Constructor
          */
         PhysicalModelBase();

         /**
          * @brief Destructor
          */
         ~PhysicalModelBase();
   };

}

#endif // PHYSICALMODELBASE_HPP
