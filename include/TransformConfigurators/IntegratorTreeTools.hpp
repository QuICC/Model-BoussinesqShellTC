/** 
 * @file IntegratorTreeTools.hpp
 * @brief Implementation of tools to work with the integrator transform trees
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef INTEGRATORTREETOOLS_HPP
#define INTEGRATORTREETOOLS_HPP

// System includes
//
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "TransformConfigurators/IntegratorBranch.hpp"
#include "TransformConfigurators/IntegratorTree.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of requirement tools for the integrator transform tree
    */
   class IntegratorTreeTools
   {
      public:
         /**
          * @brief Create vector of projector trees from map of branches
          *
          * @param rTrees     Output integrator trees
          * @param branches   Input tree branches
          */
         static void generateTrees(std::vector<IntegratorTree>& rTrees, const std::map<PhysicalNames::Id, std::vector<IntegratorBranch> >& branches);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         IntegratorTreeTools();

         /**
          * @brief Destructor
          */
         ~IntegratorTreeTools();
   };

}
}

#endif // INTEGRATORTREETOOLS_HPP
