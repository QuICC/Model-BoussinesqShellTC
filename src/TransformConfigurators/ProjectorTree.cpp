/** 
 * @file ProjectorTree.cpp
 * @brief Source of the implementation of the backward tranform tree
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
#include "TransformConfigurators/ProjectorTree.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   ProjectorTree::ProjectorTree(const PhysicalNames::Id name, const FieldComponents::Spectral::Id comp)
      :mName(name), mComp(comp)
   {
   }

   ProjectorTree::~ProjectorTree()
   {
   }

   PhysicalNames::Id ProjectorTree::name() const
   {
      return this->mName;
   }

   FieldComponents::Spectral::Id ProjectorTree::comp() const
   {
      return this->mComp;
   }

   int ProjectorTree::nSpecEdges() const
   {
      return this->mTree.size();
   }

   int ProjectorTree::nPartEdges() const
   {
      int n = 0;
      ProjectorSpecEdge_iterator it;
      for(it = this->mTree.begin(); it != this->mTree.end(); ++it)
      {
         n += it->nEdges();
      }

      return n;
   }

   ProjectorSpecEdge_range ProjectorTree::edgeRange() const
   {
      return std::make_pair(this->mTree.begin(), this->mTree.end());
   }

   ProjectorSpecEdge& ProjectorTree::addEdge(const ProjSpecId op, const int n)
   {
      this->mTree.push_back(ProjectorSpecEdge(op, n));

      return this->mTree.back();
   }

}
}
