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

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This template describes the complete integration tree for 3D space
    */
   class IntegratorTree3D
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

#endif // INTEGRATORTREE3D_HPP
