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
      :IntegratorTree2D(name, comp)
   {
   }

   IntegratorTree3D::~IntegratorTree3D()
   {
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

}
}
