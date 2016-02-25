/** 
 * @file ProjectorTree3DTools.hpp
 * @brief Implementation of tools to work with the projector transform trees for 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROJECTORTREE3DTOOLS_HPP
#define PROJECTORTREE3DTOOLS_HPP

// System includes
//
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "TransformConfigurators/ProjectorBranch3D.hpp"
#include "TransformConfigurators/ProjectorTree3D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of tools to work with the projector transform trees for 3D space
    */
   class ProjectorTree3DTools
   {
      public:
         /**
          * @brief Create vector of projector trees from map of branches
          *
          * @param rTrees     Output projector trees
          * @param branches   Input tree branches
          */
         static void generateTrees(std::vector<ProjectorTree3D>& rTrees, const std::map<PhysicalNames::Id, std::vector<ProjectorBranch3D> >& branches);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         ProjectorTree3DTools();

         /**
          * @brief Destructor
          */
         ~ProjectorTree3DTools();
   };

}
}

#endif // PROJECTORTREE3DTOOLS_HPP
