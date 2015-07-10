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
      :ProjectorTree2D(name, comp)
   {
   }

   ProjectorTree3D::~ProjectorTree3D()
   {
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

}
}
