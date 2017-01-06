/** 
 * @file TubularSplitting.hpp
 * @brief Implementation of a double load splitting algorithm, aka "Tubular" splitting
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TUBULARSPLITTING_HPP
#define TUBULARSPLITTING_HPP

// System includes
//
#include <utility>
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"
#include "LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "Resolutions/TransformResolution.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of a double load splitting algorithm, aka "Tubular" splitting
    */
   class TubularSplitting: public SplittingAlgorithm
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id ID of the CPU/Core
          * @param nCpu Number of cores used
          * @param dim  Dimensions
          */
         TubularSplitting(const int id, const int nCpu, const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~TubularSplitting();

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
          * @param status  Status output
          */
         virtual SharedTransformResolution splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status);

         /**
          * @brief Select the transform grouper
          *
          * \mhdTodo Grouper selection is not correctly implemented
          */
         void selectGrouper();

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         virtual Array computeScore(SharedResolution spResolution);

      private:
   };

}
}

#endif // TUBULARSPLITTING_HPP
