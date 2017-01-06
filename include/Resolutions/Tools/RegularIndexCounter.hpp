/** 
 * @file RegularIndexCounter.hpp
 * @brief Implementation of regular index counter 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef REGULARINDEXCOUNTER_HPP
#define REGULARINDEXCOUNTER_HPP

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
#include "Resolutions/Tools/IndexCounter.hpp"
#include "Resolutions/SimulationResolution.hpp"
#include "Resolutions/CoreResolution.hpp"

namespace QuICC {

   /**
    * @brief Implementation of regular index counter
    */ 
   class RegularIndexCounter: public IndexCounter
   {
      public:
         /**
          * @brief Constructor
          */
         RegularIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu);

         /**
          * @brief Empty destructor
          */
         ~RegularIndexCounter();

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version uses the internally stored simulation resolution
          *
          * @param spaceId Spacial the resolution represent
          */
         virtual ArrayI orderedDimensions(const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version reorders the input dimensions
          *
          * @param dims    Array of dimensions to reorder (1D, 2D, 3D, ...)
          * @param spaceId Spacial the resolution represent
          */
         virtual ArrayI orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Compute the offsets for the local modes
          */
         virtual void computeOffsets(std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Compute the offsets for the local modes by comparing to a reference simulation
          */
         virtual void computeOffsets(std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const;
         
      protected:

      private:
         /**
          * @brief Local copy of the simulation resolution
          */
         SharedCSimulationResolution mspSim;

         /**
          * @brief Local copy of the local resolution
          */
         SharedCCoreResolution mspCpu;

   };

   /// Typedef for an smart reference counting pointer for a RegularIndexCounter
   typedef SharedPtrMacro<RegularIndexCounter>   SharedRegularIndexCounter;

}

#endif // REGULARINDEXCOUNTER_HPP
