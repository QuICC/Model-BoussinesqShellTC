/** 
 * @file SerialSplitting.cpp
 * @brief Source of the implementation of a serial code "load splitting"
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

#include <iostream>
namespace GeoMHDiSCC {

namespace Parallel {

   SerialSplitting::SerialSplitting(const int id, const int nCpu, const ArrayI& dim)
      : SplittingAlgorithm(id, nCpu, dim, Splitting::Algorithms::SERIAL)
   {
      // Initialise the NCpu factors
      this->initFactors(1);
   }

   SerialSplitting::~SerialSplitting()
   {
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

   SharedTransformResolution  SerialSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId)
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
      this->mspScheme->fillIndexes(transId, fwd1D, bwd1D, idx2D, idx3D);
      std::cerr << " ------------------------- " << transId << " ----------------------" << std::endl;
      std::cerr << " ~~ 3D ~~" << std::endl;
      std::cerr << idx3D.transpose() << std::endl;
      std::cerr << " ~~ 2D ~~" << std::endl;
      for(int i = 0; i < idx3D.size(); i++)
      {
         std::cerr << idx2D.at(i).transpose() << std::endl;
      }

      // Create TransformResolution object
      return SharedTransformResolution(new TransformResolution(fwd1D, bwd1D, idx2D, idx3D));
   }

   int SerialSplitting::computeScore(SharedResolution spResolution)
   {
      return 1;
   }
}
}
