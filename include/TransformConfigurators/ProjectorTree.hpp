/** 
 * @file ProjectorTree.hpp
 * @brief This template describes the complete projection tree
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROJECTORTREE_HPP
#define PROJECTORTREE_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformSelector.hpp"
#include "TransformConfigurators/TransformEdge.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This template describes the complete projection tree
    */
   class ProjectorTree
   {  
      public:
         /// Typedefs to simplify definition of projection operators
         typedef TransformCoordinatorType::Transform1DType::ProjectorType::Id Proj1DId;
         typedef TransformCoordinatorType::Transform2DType::ProjectorType::Id Proj2DId;
         typedef TransformCoordinatorType::Transform3DType::ProjectorType::Id Proj3DId;
         typedef TransformEdge<Proj1DId,Proj2DId,Proj3DId> Projector1DEdge;
         typedef TransformEdge<Proj2DId,Proj3DId,void> Projector2DEdge;
         typedef TransformEdge<Proj3DId,void,void> Projector3DEdge;
         typedef std::vector<Projector1DEdge>::const_iterator  Projector1DEdge_iterator;
         typedef std::vector<Projector2DEdge>::const_iterator  Projector2DEdge_iterator;
         typedef std::vector<Projector3DEdge>::const_iterator  Projector3DEdge_iterator;
         typedef std::pair<Projector1DEdge_iterator,Projector1DEdge_iterator> Projector1DEdge_range;
         typedef std::pair<Projector2DEdge_iterator,Projector2DEdge_iterator> Projector2DEdge_range;
         typedef std::pair<Projector3DEdge_iterator,Projector3DEdge_iterator> Projector3DEdge_range;

         /**
          * @brief Contructor for operation
          */
         ProjectorTree(const PhysicalNames::Id name, const FieldComponents::Spectral::Id comp);

         /**
          * @brief Destructor
          */
         ~ProjectorTree();

         /**
          * @brief Get physical name
          */
         PhysicalNames::Id name() const;

         /**
          * @brief Get spectral field component
          */
         FieldComponents::Spectral::Id comp() const;

         /**
          * @brief number of 1D edges
          */
         int nEdges1D() const;

         /**
          * @brief number of 1D edges
          */
         int nEdges2D() const;

         /**
          * @brief Get vector range for 1D edges
          */
         Projector1DEdge_range edgeRange() const;

         /**
          * @brief Add an edge to tree
          */
         Projector1DEdge& addEdge(const Proj1DId op, const int n);

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
         std::vector<Projector1DEdge> mTree;
   };

}
}

#endif // PROJECTORTREE_HPP
