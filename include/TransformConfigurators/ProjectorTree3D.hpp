/** 
 * @file ProjectorTree3D.hpp
 * @brief This template describes the complete projection tree for 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROJECTORTREE3D_HPP
#define PROJECTORTREE3D_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformLeafSelector.hpp"
#include "TransformConfigurators/ProjectorTree2D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This template describes the complete projection tree for 3D space
    */
   class ProjectorTree3D: public ProjectorTree2D
   {  
      public:
         /**
          * @brief Contructor for operation
          */
         ProjectorTree3D(const PhysicalNames::Id name, const FieldComponents::Spectral::Id comp);

         /**
          * @brief Destructor
          */
         ~ProjectorTree3D();

         /**
          * @brief number of partially transformed edges
          */
         int nPartEdges() const;

      private:
   };

}
}

#endif // PROJECTORTREE3D_HPP
