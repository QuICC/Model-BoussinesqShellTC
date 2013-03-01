/** \file CommunicatorBase.hpp
 *  \brief Implementation of a communicator base
 *
 *  \mhdBug Needs test
 */

#ifndef COMMUNICATORBASE_HPP
#define COMMUNICATORBASE_HPP

// Configuration includes
//

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
    * @brief Implementation of communicator base
    */ 
   class CommunicatorBase
   {
      public:
         /**
         * @brief Constructor
         */
         CommunicatorBase();

         /**
         * @brief Destructor
         */
         virtual ~CommunicatorBase();
         
      protected:
         /**
          * @brief Create the required buffers (does no allocation yet)
          *
          * @param nBuffers Number of buffers to initialise
          */
         void createBuffers(const int nBuffers);

         /**
          * @brief Allocate buffers
          *
          * @param idx        Index of the buffer
          * @param sizes      Sizes per node
          * @param maxPacks   Maximum number of packs
          */
         void allocateBuffers(const int idx, const std::vector<int>& sizes, const int maxPacks);

         /**
          * @brief Allocate buffers
          *
          * @param idx        Index of the buffer
          * @param aSizes      Sizes per node
          * @param maxAPacks   Maximum number of packs
          * @param bSizes      Sizes per node
          * @param maxBPacks   Maximum number of packs
          */
         void allocateBuffers(const int idx, const std::vector<int>& aSizes, const int maxAPacks, const std::vector<int>& bSizes, const int maxBPacks);

         /**
          * @brief Cleanup the communication buffers
          */
         void cleanupBuffers();

         /**
          * @brief MPI communication buffers
          */
         std::vector< std::vector<char *> > mBuffers;

      private:
   };

}
}

#endif // COMMUNICATORBASE_HPP
