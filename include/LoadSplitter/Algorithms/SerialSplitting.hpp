/** 
 * @file SerialSplitting.hpp
 * @brief Implementation of a serial "load splitting" algorithm
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SERIALSPLITTING_HPP
#define SERIALSPLITTING_HPP

// System includes
//
#include <utility>
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "Resolutions/TransformResolution.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of a serial "load splitting" algorithm
    */
   class SerialSplitting: public SplittingAlgorithm
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id   ID of the CPU
          * @param nCpu Number of CPUs used
          * @param dim  Dimensions
          */
         SerialSplitting(const int id, const int nCpu, const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~SerialSplitting();

         /**
          * @brief Check if factorisation is applicable to scheme
          */
         virtual bool applicable() const;
         
      protected:
         /**
          * @brief Split ith dimension transform
          *
          * @param transId Split the ith dimension
          * @param cpuId   ID of the CPU
          */
         virtual SharedTransformResolution splitDimension(const Dimensions::Transform::Id transId, const int cpuId);

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         virtual int computeScore(SharedResolution spResolution);

      private:
   };

}
}

#endif // SERIALSPLITTING_HPP
