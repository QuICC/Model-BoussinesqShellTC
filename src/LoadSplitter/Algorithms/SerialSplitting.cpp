/** \file SerialSplitting.cpp
 *  \brief Source of the implementation of a serial code "load splitting"
 */

// System includes
//

// External includes
//

// Class include
//
#include "LoadSplitter/Algorithms/SerialSplitting.hpp"

// Project includes
//
#include "Base/IO/Ascii/FormatToolbox.hpp"

namespace GeoMHDiSCC {

   SerialSplitting::SerialSplitting(const int id, const int nCpu, const ArrayI& dim)
      : SplittingAlgorithm(id, nCpu, dim, Splittings::Algorithms::SERIAL)
   {
      // Factorise N_cpu
      this->factoriseNCpu(1);

      // Filter factors
      this->filterFactors();
   }

   bool SerialSplitting::applicable() const
   {
      bool status = true;

      // As long as there is a single CPU it should be applicable
      status = (status && (this->nCpu() == 1));

      // Check for scheme specific conditions
      status = (status && this->mspScheme->applicable());

      return status;
   }

   SharedTransformResolution  SerialSplitting::splitDimension(const int dim, const int id)
   {
      // Storage for the forward 1D indexes
      std::vector<ArrayI>  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<ArrayI>  bwd1D;
      // Storage for the 2D indexes
      std::vector<ArrayI>  idx2D;
      // Storage for the 3D indexes
      ArrayI  idx3D;

      // Compute the indexes
      this->mspScheme->fillIndexes(dim, fwd1D, bwd1D, idx2D, idx3D);

      // Create TransformResolution object
      return SharedTransformResolution(new TransformResolution(fwd1D, bwd1D, idx2D, idx3D));
   }

   int SerialSplitting::computeScore(SharedResolution spResolution)
   {
      return 1;
   }
}
