/** \file SpatialScheme.cpp
 *  \brief Source of the base for the scheme implementations
 */

// System includes
//

// External includes
//

// Class include
//
#include "Base/SpatialSchemes/SpatialScheme.hpp"

// Project includes
//

namespace EPMPhoenix {

   SpatialScheme::SpatialScheme(const int dims)
      : SchemeBase(dims), mDims(dims)
   {
   }

   SpatialScheme::~SpatialScheme()
   {
   }

   void SpatialScheme::init()
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
      this->setMemory();
   }

}
