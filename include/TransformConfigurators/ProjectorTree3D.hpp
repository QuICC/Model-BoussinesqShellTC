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

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This template describes the complete projection tree for 3D space
    */
   class ProjectorTree3D
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
          * @brief Get physical name
          */
         PhysicalNames::Id name() const;

         /**
          * @brief Get spectral field component
          */
         FieldComponents::Spectral::Id comp() const;

         /**
          * @brief number of spectral space edges
          */
         int nSpecEdges() const;

         /**
          * @brief number of partially transformed edges
          */
         int nPartEdges() const;

         /**
          * @brief Get vector range for spectral edges
          */
         ProjectorSpecEdge_range edgeRange() const;

         /**
          * @brief Add a spectral edge to tree
          */
         ProjectorSpecEdge& addEdge(const ProjSpecId op, const int n);

      private:
         /**
          * @brief Physical field name
          */
         PhysicalNames::Id mName;

         /**
          * @brief Spectral field component
          */
         FieldComponents::Spectral::Id mComp;

         /**
          * @brief Projector edges of the tree
          */
         std::vector<ProjectorSpecEdge> mTree;
   };

}
}

#endif // PROJECTORTREE3D_HPP
