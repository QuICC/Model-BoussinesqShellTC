/** 
 * @file SimulationBoundary.hpp
 * @brief Implementation of a simple simulation wide boundary condition interface
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SIMULATIONBOUNDARY_HPP
#define SIMULATIONBOUNDARY_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "SmartPointers/SharedPtrMacro.h"
#include "Enums/FieldIds.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a simple simulation wide boundary condition interface
    */
   class SimulationBoundary
   {
      public:
         /**
          * @brief Constructor
          */
         SimulationBoundary(const std::map<std::string,int>& bcIds);

         /**
          * @brief Destructor
          */
         ~SimulationBoundary();

         /**
          * @brief Get tag map
          */
         std::map<std::string,int> getTagMap() const;

         /**
          * @brief Get tag map
          */
         int bcId(const PhysicalNames::Id id) const;
         
      protected:

      private:
         /**
          * @brief Convert tag map to ID map
          *
          * @param bcIds   Tag map
          */
         void convert(const std::map<std::string,int>& bcIds);

         /**
          * @brief Storage for the boundary conditions
          */
         std::map<PhysicalNames::Id, int> mBcs;
   };

   /// Typedef for a shared pointer to a SimulationBoundary object
   typedef SharedPtrMacro<SimulationBoundary>   SharedSimulationBoundary;
}

#endif // SIMULATIONBOUNDARY_HPP
