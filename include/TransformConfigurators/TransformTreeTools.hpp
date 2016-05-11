/** 
 * @file TransformTreeTools.hpp
 * @brief Implementation of tools to work with the transform trees
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMTREETOOLS_HPP
#define TRANSFORMTREETOOLS_HPP

// System includes
//
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "TransformConfigurators/TransformPath.hpp"
#include "TransformConfigurators/TransformTree.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Implementation of tools to work with the projector transform trees for 3D space
    */
   class TransformTreeTools
   {
      public:
         /**
          * @brief Create vector of transform trees from vector of transform paths
          *
          * @param rTrees     Output transform trees
          * @param branches   Input transfrom paths
          */
         static void generateTrees(std::vector<TransformTree>& rTrees, const std::map<PhysicalNames::Id, std::vector<TransformPath> >& branches);
         
      protected:

      private:
         /**
          * @brief Recursive grows of the tree
          */
         static void growTree(std::map<PhysicalNames::Id, std::vector<TransformPath> >::const_iterator nameIt, std::set<int>::const_iterator compIt, const int dim, std::vector<int>& path, TransformTreeEdge& rPrev);

         /**
          * @brief Constructor
          */
         TransformTreeTools();

         /**
          * @brief Destructor
          */
         ~TransformTreeTools();
   };

}
}

#endif // TRANSFORMTREETOOLS_HPP
