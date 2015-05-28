/** 
 * @file ProjectorTreeTools.hpp
 * @brief Implementation of tools to work with the projector transform trees
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROJECTORTREETOOLS_HPP
#define PROJECTORTREETOOLS_HPP

// System includes
//
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "TransformConfigurators/ProjectorBranch.hpp"
#include "TransformConfigurators/ProjectorTree.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of requirement tools for equations and variables
    */
   class ProjectorTreeTools
   {
      public:
         /**
          * @brief Create vector of projector trees from map of branches
          *
          * @param rTrees     Output projector trees
          * @param branches   Input tree branches
          */
         static void generateTrees(std::vector<ProjectorTree>& rTrees, const std::map<PhysicalNames::Id, std::vector<ProjectorBranch> >& branches);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         ProjectorTreeTools();

         /**
          * @brief Destructor
          */
         ~ProjectorTreeTools();
   };

}
}

#endif // PROJECTORTREETOOLS_HPP
