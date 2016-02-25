/** 
 * @file ProjectorTree2D.hpp
 * @brief This template describes the complete projection tree for 2D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROJECTORTREE2D_HPP
#define PROJECTORTREE2D_HPP

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
    * @brief This template describes the complete projection tree for 2D space
    */
   class ProjectorTree2D
   {  
      public:
         /**
          * @brief Contructor for operation
          */
         ProjectorTree2D(const PhysicalNames::Id name, const FieldComponents::Spectral::Id comp);

         /**
          * @brief Destructor
          */
         ~ProjectorTree2D();

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
          * @brief Get vector range for spectral edges
          */
         ProjectorSpecEdge_range edgeRange() const;

         /**
          * @brief Add a spectral edge to tree
          */
         ProjectorSpecEdge& addEdge(const ProjSpecId op, const int n);

      protected:
         /**
          * @brief Projector edges of the tree
          */
         std::vector<ProjectorSpecEdge> mTree;

      private:
         /**
          * @brief Physical field name
          */
         PhysicalNames::Id mName;

         /**
          * @brief Spectral field component
          */
         FieldComponents::Spectral::Id mComp;
   };

}
}

#endif // PROJECTORTREE2D_HPP
