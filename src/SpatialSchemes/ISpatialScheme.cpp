/** 
 * @file ISpatialScheme.cpp
 * @brief Source of the base for the scheme implementations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <vector>

// External includes
//

// Class include
//
#include "SpatialSchemes/ISpatialScheme.hpp"

// Project includes
//
#include "Resolutions/Tools/RegularIndexCounter.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   void ISpatialScheme::interpretConfigDimensions(ArrayI& rDim)
   {
   }

   void ISpatialScheme::tuneResolution(SharedResolution spRes)
   {
      SharedRegularIndexCounter   spCounter(new RegularIndexCounter(spRes->sim(), spRes->cpu()));

      spRes->setIndexCounter(spCounter);

      ISpatialScheme::tuneMpiResolution();
   }

   ISpatialScheme::ISpatialScheme(const int dims)
      : ISchemeCosts(dims), mDims(dims)
   {
   }

   ISpatialScheme::~ISpatialScheme()
   {
   }

   void ISpatialScheme::init()
   {
      // Initialise Storage for the dimensions
      for(int i = 0; i < this->dims(); i++)
      {
         this->mDimensions.push_back(ArrayI(this->dims()+1));
      }

      // Initialise the domain's dimensions
      this->setDimensions(); 

      // Set the transform costs
      this->setCosts();

      // Set the transform scalings
      this->setScalings();

      // Set the memory costs
      this->setMemoryScore();
   }

   int ISpatialScheme::dim(const Dimensions::Transform::Id transId, const Dimensions::Data::Id dataId) const
   {
      // Assert for domain dimensions
      assert(static_cast<int>(transId) < this->mDims);

      // Assert for dimension size
      assert(this->mDimensions.at(static_cast<int>(transId)).size() > static_cast<int>(dataId));

      return this->mDimensions.at(static_cast<int>(transId))(static_cast<int>(dataId));
   }

   void ISpatialScheme::setTransformSpace(const ArrayI& dim)
   {
      this->mTransformSpace = dim;
   }

   const ArrayI& ISpatialScheme::getTransformSpace() const
   {
      assert(this->mTransformSpace.size() == this->dims());
      assert(this->mTransformSpace.array().abs().minCoeff() > 0);

      return this->mTransformSpace;
   }

   void ISpatialScheme::setDimension(int n, const Dimensions::Transform::Id transId, const Dimensions::Data::Id dataId)
   {
      // Assert for positive size
      assert(n > 0);

      // Assert for domain dimensions
      assert(static_cast<int>(transId) < this->mDims);

      // Assert for dimension size
      assert(this->mDimensions.at(static_cast<int>(transId)).size() > static_cast<int>(dataId));

      this->mDimensions.at(static_cast<int>(transId))(static_cast<int>(dataId)) = n;
   }

   int ISpatialScheme::dims() const
   {
      return this->mDims;
   }

   void ISpatialScheme::tuneMpiResolution()
   {
      // Create single rank communicator
      #ifdef GEOMHDISCC_MPI
         std::vector<int>  ranks;
         ranks.push_back(FrameworkMacro::id());
         FrameworkMacro::setSpectralComm(ranks);
      #endif //GEOMHDISCC_MPI
   }

}
}
