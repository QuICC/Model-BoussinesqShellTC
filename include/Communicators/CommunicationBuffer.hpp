/**
 * @file CommunicationBuffer.hpp
 * @brief Implementation of a "raw" communication buffer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef COMMUNICATIONBUFFER_HPP
#define COMMUNICATIONBUFFER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * 
    * @brief Implementation of a "raw" communication buffer
    */ 
   class CommunicationBuffer
   {
      public:
         /**
         * @brief Constructor
         */
         CommunicationBuffer();

         /**
         * @brief Destructor
         */
         ~CommunicationBuffer();

         /**
          * @brief Allocate buffers
          *
          * @param sizes      Sizes per node
          * @param maxPacks   Maximum number of packs
          */
         void allocate(const std::vector<int>& sizes, const int maxPacks);

         /**
          * @brief Allocate buffers to fit both requested sizes
          *
          * @param aSizes      Sizes per node
          * @param maxAPacks   Maximum number of packs
          * @param bSizes      Sizes per node
          * @param maxBPacks   Maximum number of packs
          */
         void allocateMax(const std::vector<int>& aSizes, const int maxAPacks, const std::vector<int>& bSizes, const int maxBPacks);

         /**
          * @brief Get pointer to raw buffer storage
          */
         char* at(const int id);
         
      protected:

      private:
         /**
          * @brief Cleanup the communication buffers
          */
         void cleanupBuffers();

         /**
          * @brief MPI communication buffers
          */
         std::vector<char *> mData;
   };

   /// Typedef for a shared pointer to a converter buffer object
   typedef SharedPtrMacro<CommunicationBuffer>   SharedCommunicationBuffer;

}
}

#endif // COMMUNICATIONBUFFER_HPP
