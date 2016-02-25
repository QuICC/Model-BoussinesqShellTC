/** 
 * @file ProjectorTree2D.cpp
 * @brief Source of the implementation of the backward tranform tree for 2D space
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
#include "TransformConfigurators/ProjectorTree2D.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   ProjectorTree2D::ProjectorTree2D(const PhysicalNames::Id name, const FieldComponents::Spectral::Id comp)
      :mName(name), mComp(comp)
   {
   }

   ProjectorTree2D::~ProjectorTree2D()
   {
   }

   PhysicalNames::Id ProjectorTree2D::name() const
   {
      return this->mName;
   }

   FieldComponents::Spectral::Id ProjectorTree2D::comp() const
   {
      return this->mComp;
   }

   int ProjectorTree2D::nSpecEdges() const
   {
      return this->mTree.size();
   }

   ProjectorSpecEdge_range ProjectorTree2D::edgeRange() const
   {
      return std::make_pair(this->mTree.begin(), this->mTree.end());
   }

   ProjectorSpecEdge& ProjectorTree2D::addEdge(const ProjSpecId op, const int n)
   {
      this->mTree.push_back(ProjectorSpecEdge(op, n));

      return this->mTree.back();
   }

}
}
