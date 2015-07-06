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
#include "TypeSelectors/TreeSelector.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This template describes the complete integration tree
    */
   class IntegratorTree
   {  
      public:
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
          * @brief number of physical space edges
          */
         int nPhysEdges() const;

         /**
          * @brief number of partially transformed edges
          */
         int nPartEdges() const;

         /**
          * @brief Get vector range for physical space edges
          */
         IntegratorPhysEdge_range edgeRange() const;

         /**
          * @brief Add an edge to tree
          */
         IntegratorPhysEdge& addEdge(const IntgPhysId op, const int n);

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
         std::vector<IntegratorPhysEdge> mTree;
   };

}
}

#endif // INTEGRATORTREE_HPP
