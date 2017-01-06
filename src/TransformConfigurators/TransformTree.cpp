/** 
 * @file TransformTree.cpp
 * @brief Source of the implementation of the tranform tree
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/TransformTree.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

   TransformTree::TransformTree(const PhysicalNames::Id name, const int comp)
      :mName(name), mComp(comp), mRoot(-1, 1)
   {
   }

   TransformTree::~TransformTree()
   {
   }

   PhysicalNames::Id TransformTree::name() const
   {
      return this->mName;
   }

   const TransformTreeEdge& TransformTree::root() const
   {
      return this->mRoot;
   }

   TransformTreeEdge& TransformTree::rRoot()
   {
      return this->mRoot;
   }

   int TransformTree::nEdges(const int depth) const
   {
      return this->mRoot.nEdges(depth);
   }

}
}
