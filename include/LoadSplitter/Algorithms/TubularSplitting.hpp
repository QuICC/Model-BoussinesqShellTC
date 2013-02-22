/** \file TubularSplitting.hpp
 *  \brief Implementation of a double load splitting algorithm, aka "Tubular" splitting
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
#include "LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "Base/Resolutions/TransformResolution.hpp"

namespace GeoMHDiSCC {

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
          * @param dim  Split the ith transform
          * @param id   ID of the CPU
          */
         virtual SharedTransformResolution splitDimension(const int dim, const int id);

         /**
          * @brief Select the transform grouper
          */
         void selectGrouper();

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         virtual int computeScore(SharedResolution spResolution);

      private:
   };

}

#endif // TUBULARSPLITTING_HPP
