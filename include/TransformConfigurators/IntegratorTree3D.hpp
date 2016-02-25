/** 
 * @file IntegratorTree3D.hpp
 * @brief This template describes the complete integration tree for 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef INTEGRATORTREE3D_HPP
#define INTEGRATORTREE3D_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformLeafSelector.hpp"
#include "TransformConfigurators/IntegratorTree2D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This template describes the complete integration tree for 3D space
    */
   class IntegratorTree3D: public IntegratorTree2D
   {  
      public:
         /**
          * @brief Contructor for operation
          */
         IntegratorTree3D(const PhysicalNames::Id name, const FieldComponents::Physical::Id comp);

         /**
          * @brief Destructor
          */
         ~IntegratorTree3D();

         /**
          * @brief number of partially transformed edges
          */
         int nPartEdges() const;

      private:
   };

}
}

#endif // INTEGRATORTREE3D_HPP
