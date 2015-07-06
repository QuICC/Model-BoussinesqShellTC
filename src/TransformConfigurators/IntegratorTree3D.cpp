/** 
 * @file IntegratorTree3D.cpp
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
#include "TransformConfigurators/IntegratorTree3D.hpp"

// Intgect includes
//

namespace GeoMHDiSCC {

namespace Transform {

   IntegratorTree3D::IntegratorTree3D(const PhysicalNames::Id name, const FieldComponents::Physical::Id comp)
      :mName(name), mComp(comp)
   {
   }

   IntegratorTree3D::~IntegratorTree3D()
   {
   }

   PhysicalNames::Id IntegratorTree3D::name() const
   {
      return this->mName;
   }

   FieldComponents::Physical::Id IntegratorTree3D::comp() const
   {
      return this->mComp;
   }

   int IntegratorTree3D::nPhysEdges() const
   {
      return this->mTree.size();
   }

   int IntegratorTree3D::nPartEdges() const
   {
      int n = 0;
      IntegratorPhysEdge_iterator it;
      for(it = this->mTree.begin(); it != this->mTree.end(); ++it)
      {
         n += it->nEdges();
      }

      return n;
   }

   IntegratorPhysEdge_range IntegratorTree3D::edgeRange() const
   {
      return std::make_pair(this->mTree.begin(), this->mTree.end());
   }

   IntegratorPhysEdge& IntegratorTree3D::addEdge(const IntgPhysId op, const int n)
   {
      this->mTree.push_back(IntegratorPhysEdge(op, n));

      return this->mTree.back();
   }

}
}
