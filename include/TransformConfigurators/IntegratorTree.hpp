/** 
 * @file IntegratorTree.hpp
 * @brief This template describes the complete integration tree
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef INTEGRATORTREE_HPP
#define INTEGRATORTREE_HPP

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
    * @brief This template describes the complete integration tree
    */
   class IntegratorTree
   {  
      public:
         /// Typedefs to simplify definition of integration operators
         typedef TransformCoordinatorType::Transform1DType::IntegratorType::Id Intg1DId;
         typedef TransformCoordinatorType::Transform2DType::IntegratorType::Id Intg2DId;
         typedef TransformCoordinatorType::Transform3DType::IntegratorType::Id Intg3DId;
         typedef TransformEdge<Intg3DId,Intg2DId,Intg1DId> Integrator3DEdge;
         typedef TransformEdge<Intg2DId,Intg1DId,void> Integrator2DEdge;
         typedef TransformEdge<Intg1DId,void,void> Integrator1DEdge;
         typedef std::vector<Integrator1DEdge>::const_iterator  Integrator1DEdge_iterator;
         typedef std::vector<Integrator2DEdge>::const_iterator  Integrator2DEdge_iterator;
         typedef std::vector<Integrator3DEdge>::const_iterator  Integrator3DEdge_iterator;
         typedef std::pair<Integrator1DEdge_iterator,Integrator1DEdge_iterator> Integrator1DEdge_range;
         typedef std::pair<Integrator2DEdge_iterator,Integrator2DEdge_iterator> Integrator2DEdge_range;
         typedef std::pair<Integrator3DEdge_iterator,Integrator3DEdge_iterator> Integrator3DEdge_range;

         /**
          * @brief Contructor for operation
          */
         IntegratorTree(const PhysicalNames::Id name, const FieldComponents::Physical::Id comp);

         /**
          * @brief Destructor
          */
         ~IntegratorTree();

         /**
          * @brief Get physical name
          */
         PhysicalNames::Id name() const;

         /**
          * @brief Get spectral field component
          */
         FieldComponents::Physical::Id comp() const;

         /**
          * @brief number of 3D edges
          */
         int nEdges3D() const;

         /**
          * @brief number of 1D edges
          */
         int nEdges2D() const;

         /**
          * @brief Get vector range for 3D edges
          */
         Integrator3DEdge_range edgeRange() const;

         /**
          * @brief Add an edge to tree
          */
         Integrator3DEdge& addEdge(const Intg3DId op, const int n);

      private:
         /**
          * @brief Physical field name
          */
         PhysicalNames::Id mName;

         /**
          * @brief Physical field component
          */
         FieldComponents::Physical::Id mComp;

         /**
          * @brief Integrator edges of the tree
          */
         std::vector<Integrator3DEdge> mTree;
   };

}
}

#endif // INTEGRATORTREE_HPP
