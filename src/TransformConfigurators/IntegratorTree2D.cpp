/** 
 * @file IntegratorTree2D.cpp
 * @brief Source of the implementation of the forward transform tree for 3D space
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
#include "TransformConfigurators/IntegratorTree2D.hpp"

// Intgect includes
//

namespace GeoMHDiSCC {

namespace Transform {

   IntegratorTree2D::IntegratorTree2D(const PhysicalNames::Id name, const FieldComponents::Physical::Id comp)
      :mName(name), mComp(comp)
   {
   }

   IntegratorTree2D::~IntegratorTree2D()
   {
   }

   PhysicalNames::Id IntegratorTree2D::name() const
   {
      return this->mName;
   }

   FieldComponents::Physical::Id IntegratorTree2D::comp() const
   {
      return this->mComp;
   }

   int IntegratorTree2D::nPhysEdges() const
   {
      return this->mTree.size();
   }

   IntegratorPhysEdge_range IntegratorTree2D::edgeRange() const
   {
      return std::make_pair(this->mTree.begin(), this->mTree.end());
   }

   IntegratorPhysEdge& IntegratorTree2D::addEdge(const IntgPhysId op, const int n)
   {
      this->mTree.push_back(IntegratorPhysEdge(op, n));

      return this->mTree.back();
   }

}
}
