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

   void ISpatialScheme::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      ISpatialScheme::tuneMpiResolution(descr);
   }

   ISpatialScheme::ISpatialScheme(const int dims)
      : ISchemeCosts(dims), mDims(dims)
   {
   }

   ISpatialScheme::~ISpatialScheme()
   {
   }

   void ISpatialScheme::addIndexCounter(SharedResolution spRes)
   {
      SharedRegularIndexCounter   spCounter(new RegularIndexCounter(spRes->sim(), spRes->cpu()));

      spRes->setIndexCounter(spCounter);

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

   void ISpatialScheme::tuneMpiResolution(const Parallel::SplittingDescription& descr)
   {
      #if defined QUICC_MPI
         FrameworkMacro::initTransformComm(descr.structure.size());
         for(std::vector<std::multimap<int,int> >::const_iterator vIt = descr.structure.begin(); vIt != descr.structure.end(); ++vIt)
         {
            // Extract the communication group from structure
            std::multimap<int,int>::const_iterator it;
            std::set<int> filter;
            filter.insert(FrameworkMacro::id());
            for(it = vIt->equal_range(FrameworkMacro::id()).first; it != vIt->equal_range(FrameworkMacro::id()).second; ++it)
            {
               filter.insert(it->second);
            }

            // Convert set to array of CPUs in group
            ArrayI groupCpu(filter.size());
            int i = 0;
            for(std::set<int>::iterator it = filter.begin(); it != filter.end(); ++it)
            {
               groupCpu(i) = *it;
               ++i;
            }

            FrameworkMacro::addTransformComm(groupCpu);

            // Synchronize
            FrameworkMacro::synchronize();
         }

      #endif //defined QUICC_MPI
   }

}
}
