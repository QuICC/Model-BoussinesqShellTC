/** \file CommunicationBuffer.hpp
 *  \brief Implementation of a buffer for the MPI converter
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
    * @brief Implementation of a bffer for the MPI converter
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
         virtual ~CommunicationBuffer();

         /**
          * @brief Allocate buffers
          *
          * @param sizes      Sizes per node
          * @param maxPacks   Maximum number of packs
          */
         void allocate(const std::vector<int>& sizes, const int maxPacks);

         /**
          * @brief Allocate buffers
          *
          * @param aSizes      Sizes per node
          * @param maxAPacks   Maximum number of packs
          * @param bSizes      Sizes per node
          * @param maxBPacks   Maximum number of packs
          */
         void allocate(const std::vector<int>& aSizes, const int maxAPacks, const std::vector<int>& bSizes, const int maxBPacks);
         
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
