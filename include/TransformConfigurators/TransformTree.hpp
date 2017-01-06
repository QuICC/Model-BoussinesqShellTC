/** 
 * @file TransformTree.hpp
 * @brief This template describes the complete projection tree for 2D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMTREE_HPP
#define TRANSFORMTREE_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TransformConfigurators/TransformTreeEdge.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This template describes the complete projection tree for 2D space
    */
   class TransformTree
   {  
      public:
         /**
          * @brief Contructor for operation
          */
         TransformTree(const PhysicalNames::Id name, const int comp);

         /**
          * @brief Destructor
          */
         ~TransformTree();

         /**
          * @brief Get physical name
          */
         PhysicalNames::Id name() const;

         /**
          * @brief Get field component
          */
         template <typename TId> TId comp() const;

         /**
          * @brief Number of edges at given depth
          */
         int nEdges(const int depth) const;

         /**
          * @brief Get root of the tree
          */
         const TransformTreeEdge& root() const;

         /**
          * @brief Set root of the tree
          */
         TransformTreeEdge& rRoot();

      protected:

      private:
         /**
          * @brief Physical field name
          */
         PhysicalNames::Id mName;

         /**
          * @brief Root field component
          */
         int mComp;

         /**
          * @brief Root of the tree
          */
         TransformTreeEdge mRoot;
   };

   template <typename TId> inline TId TransformTree::comp() const
   {
      return static_cast<TId>(this->mComp);
   }

}
}

#endif // TRANSFORMTREE_HPP
