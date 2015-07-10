/** 
 * @file TransformTreeSelector.hpp
 * @brief Typedefs to setup the transform trees
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMTREESELECTOR_HPP
#define TRANSFORMTREESELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#ifdef GEOMHDISCC_SPATIALDIMENSION_3D
   #include "TransformConfigurators/ProjectorBranch3D.hpp"
   #include "TransformConfigurators/ProjectorTree3D.hpp"
   #include "TransformConfigurators/ProjectorTree3DTools.hpp"

   #include "TransformConfigurators/IntegratorBranch3D.hpp"
   #include "TransformConfigurators/IntegratorTree3D.hpp"
   #include "TransformConfigurators/IntegratorTree3DTools.hpp"
#else
   #include "TransformConfigurators/ProjectorBranch2D.hpp"
   #include "TransformConfigurators/ProjectorTree2D.hpp"
   #include "TransformConfigurators/ProjectorTree2DTools.hpp"

   #include "TransformConfigurators/IntegratorBranch2D.hpp"
   #include "TransformConfigurators/IntegratorTree2D.hpp"
   #include "TransformConfigurators/IntegratorTree2DTools.hpp"
#endif //GEOMHDISCC_SPATIALDIMENSION_3D

namespace GeoMHDiSCC {

   namespace Transform {

      #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
         typedef ProjectorBranch3D  ProjectorBranch;
         typedef ProjectorTree3D  ProjectorTree;
         typedef ProjectorTree3DTools  ProjectorTreeTools;
         
         typedef IntegratorBranch3D  IntegratorBranch;
         typedef IntegratorTree3D  IntegratorTree;
         typedef IntegratorTree3DTools  IntegratorTreeTools;
      #else
         typedef ProjectorBranch2D  ProjectorBranch;
         typedef ProjectorTree2D  ProjectorTree;
         typedef ProjectorTree2DTools  ProjectorTreeTools;
         
         typedef IntegratorBranch2D  IntegratorBranch;
         typedef IntegratorTree2D  IntegratorTree;
         typedef IntegratorTree2DTools  IntegratorTreeTools;
      #endif //GEOMHDISCC_SPATIALDIMENSION_3D
   }
}

#endif // TRANSFORMTREESELECTOR_HPP
