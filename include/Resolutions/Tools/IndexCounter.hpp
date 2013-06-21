/** \file IndexCounter.hpp
 *  \brief Implementation of base class for a generalized index counter
 */

#ifndef INDEXCOUNTER_HPP
#define INDEXCOUNTER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"

namespace GeoMHDiSCC {

   // Forward declaration to avoid header loop
   class SimulationResolution;
   typedef SharedPtrMacro<const SimulationResolution> SharedCSimulationResolution;

   /**
    * @brief Implementation of base class for a generalized index counter
    */ 
   class IndexCounter
   {
      public:
         /// Typedef for a very large unsigned int for the offsets (to avoid including hdf5.h)
         typedef long long unsigned int OffsetType;
         /**
          * @brief Constructor
          */
         IndexCounter();

         /**
          * @brief Empty destructor
          */
         ~IndexCounter();

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version uses the internally stored simulation resolution
          *
          * @param spaceId Spacial the resolution represent
          */
         virtual ArrayI orderedDimensions(const Dimensions::Space::Id spaceId) const = 0;

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version reorders the input dimensions
          *
          * @param dims    Array of dimensions to reorder (1D, 2D, 3D, ...)
          * @param spaceId Spacial the resolution represent
          */
         virtual ArrayI orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const = 0;

         /**
          * @brief Compute the offset for local modes
          */
         virtual void computeOffsets(std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const = 0;

         /**
          * @brief Compute the offset for local modes by comparing to a reference simulation
          */
         virtual void computeOffsets(std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const = 0;
         
      protected:

      private:

   };

   /// Typedef for an smart reference counting pointer for a IndexCounter
   typedef SharedPtrMacro<IndexCounter>   SharedIndexCounter;

}

#endif // INDEXCOUNTER_HPP
