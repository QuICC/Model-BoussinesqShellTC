/** 
 * @file IntegratorTree.cpp
 * @brief Source of the implementation of the forward transform tree
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
#include "TransformConfigurators/IntegratorTree.hpp"

// Intgect includes
//

namespace GeoMHDiSCC {

namespace Transform {

   IntegratorTree::IntegratorTree(const PhysicalNames::Id name, const FieldComponents::Physical::Id comp)
      :mName(name), mComp(comp)
   {
   }

   IntegratorTree::~IntegratorTree()
   {
   }

   PhysicalNames::Id IntegratorTree::name() const
   {
      return this->mName;
   }

   FieldComponents::Physical::Id IntegratorTree::comp() const
   {
      return this->mComp;
   }

   int IntegratorTree::nPhysEdges() const
   {
      return this->mTree.size();
   }

   int IntegratorTree::nPartEdges() const
   {
      int n = 0;
      IntegratorPhysEdge_iterator it;
      for(it = this->mTree.begin(); it != this->mTree.end(); ++it)
      {
         n += it->nEdges();
      }

      return n;
   }

   IntegratorPhysEdge_range IntegratorTree::edgeRange() const
   {
      return std::make_pair(this->mTree.begin(), this->mTree.end());
   }

   IntegratorPhysEdge& IntegratorTree::addEdge(const IntgPhysId op, const int n)
   {
      this->mTree.push_back(IntegratorPhysEdge(op, n));

      return this->mTree.back();
   }

}
}
