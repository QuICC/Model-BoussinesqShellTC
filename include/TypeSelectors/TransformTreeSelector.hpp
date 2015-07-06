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
#if defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TT || defined GEOMHDISCC_SPATIALSCHEME_AF || defined GEOMHDISCC_SPATIALSCHEME_CF
   #include "TransformConfigurators/ProjectorBranch2D.hpp"
   #include "TransformConfigurators/ProjectorTree2D.hpp"
   #include "TransformConfigurators/ProjectorTree2DTools.hpp"

   #include "TransformConfigurators/IntegratorBranch2D.hpp"
   #include "TransformConfigurators/IntegratorTree2D.hpp"
   #include "TransformConfigurators/IntegratorTree2DTools.hpp"
#else
   #include "TransformConfigurators/ProjectorBranch3D.hpp"
   #include "TransformConfigurators/ProjectorTree3D.hpp"
   #include "TransformConfigurators/ProjectorTree3DTools.hpp"

   #include "TransformConfigurators/IntegratorBranch3D.hpp"
   #include "TransformConfigurators/IntegratorTree3D.hpp"
   #include "TransformConfigurators/IntegratorTree3DTools.hpp"
#endif //defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TT || defined GEOMHDISCC_SPATIALSCHEME_AF || defined GEOMHDISCC_SPATIALSCHEME_CF

namespace GeoMHDiSCC {

   namespace Transform {

      #if defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TT || defined GEOMHDISCC_SPATIALSCHEME_AF || defined GEOMHDISCC_SPATIALSCHEME_CF
         typedef ProjectorBranch2D  ProjectorBranch;
         typedef ProjectorTree2D  ProjectorTree;
         typedef ProjectorTree2DTools  ProjectorTreeTools;
         
         typedef IntegratorBranch2D  IntegratorBranch;
         typedef IntegratorTree2D  IntegratorTree;
         typedef IntegratorTree2DTools  IntegratorTreeTools;
      #else
         typedef ProjectorBranch3D  ProjectorBranch;
         typedef ProjectorTree3D  ProjectorTree;
         typedef ProjectorTree3DTools  ProjectorTreeTools;
         
         typedef IntegratorBranch3D  IntegratorBranch;
         typedef IntegratorTree3D  IntegratorTree;
         typedef IntegratorTree3DTools  IntegratorTreeTools;
      #endif //defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TT || defined GEOMHDISCC_SPATIALSCHEME_AF || defined GEOMHDISCC_SPATIALSCHEME_CF
   }
}

#endif // TRANSFORMTREESELECTOR_HPP
