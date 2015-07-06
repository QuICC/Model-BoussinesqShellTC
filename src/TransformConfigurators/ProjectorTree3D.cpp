/** 
 * @file ProjectorTree3D.cpp
 * @brief Source of the implementation of the backward tranform tree for 3D space
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
#include "TransformConfigurators/ProjectorTree3D.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   ProjectorTree3D::ProjectorTree3D(const PhysicalNames::Id name, const FieldComponents::Spectral::Id comp)
      :mName(name), mComp(comp)
   {
   }

   ProjectorTree3D::~ProjectorTree3D()
   {
   }

   PhysicalNames::Id ProjectorTree3D::name() const
   {
      return this->mName;
   }

   FieldComponents::Spectral::Id ProjectorTree3D::comp() const
   {
      return this->mComp;
   }

   int ProjectorTree3D::nSpecEdges() const
   {
      return this->mTree.size();
   }

   int ProjectorTree3D::nPartEdges() const
   {
      int n = 0;
      ProjectorSpecEdge_iterator it;
      for(it = this->mTree.begin(); it != this->mTree.end(); ++it)
      {
         n += it->nEdges();
      }

      return n;
   }

   ProjectorSpecEdge_range ProjectorTree3D::edgeRange() const
   {
      return std::make_pair(this->mTree.begin(), this->mTree.end());
   }

   ProjectorSpecEdge& ProjectorTree3D::addEdge(const ProjSpecId op, const int n)
   {
      this->mTree.push_back(ProjectorSpecEdge(op, n));

      return this->mTree.back();
   }

}
}
