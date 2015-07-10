/** 
 * @file ProjectorTree2DTools.hpp
 * @brief Implementation of tools to work with the projector transform trees for 2D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROJECTORTREE2DTOOLS_HPP
#define PROJECTORTREE2DTOOLS_HPP

// System includes
//
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "TransformConfigurators/ProjectorBranch2D.hpp"
#include "TransformConfigurators/ProjectorTree2D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of tools to work with the projector transform trees for 2D space
    */
   class ProjectorTree2DTools
   {
      public:
         /**
          * @brief Create vector of projector trees from map of branches
          *
          * @param rTrees     Output projector trees
          * @param branches   Input tree branches
          */
         static void generateTrees(std::vector<ProjectorTree2D>& rTrees, const std::map<PhysicalNames::Id, std::vector<ProjectorBranch2D> >& branches);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         ProjectorTree2DTools();

         /**
          * @brief Destructor
          */
         ~ProjectorTree2DTools();
   };

}
}

#endif // PROJECTORTREE2DTOOLS_HPP
