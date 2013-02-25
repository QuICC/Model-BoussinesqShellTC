/** \file SerialSplitting.hpp
 *  \brief Implementation of a serial "load splitting" algorithm
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

   /**
    * \brief Implementation of a serial "load splitting" algorithm
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
          * @param dim  Split the ith transform
          * @param id   ID of the CPU
          */
         virtual SharedTransformResolution splitDimension(const int dim, const int id);

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         virtual int computeScore(SharedResolution spResolution);

      private:
   };

}

#endif // SERIALSPLITTING_HPP
