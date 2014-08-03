/** 
 * @file FixedSplitting.cpp
 * @brief Source of the implementation of a fixed load splitting algorithm
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <utility>
#include <queue>
#include <set>
#include <map>

// External includes
//

// Class include
//
#include "LoadSplitter/Algorithms/FixedSplitting.hpp"

// Project includes
//
#include "LoadSplitter/Algorithms/SplittingTools.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   FixedSplitting::FixedSplitting(const int id, const int nCpu, const ArrayI& dim)
      : SplittingAlgorithm(id, nCpu, dim, Splitting::Algorithms::FIXED)
   {
      // Initialise the NCpu factors
      this->initFactors(1);
   }

   FixedSplitting::~FixedSplitting()
   {
   }

   bool FixedSplitting::applicable() const
   {
      bool status = true;

      // Only works with at least 2 dimensions
      status = (status && (this->dims() > 1));

      // Check for scheme specific conditions
      status = (status && this->mspScheme->applicable());

      return status;
   }

   SharedTransformResolution  FixedSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId)
   {
      // Get size of the splittable dimension(s)
      int tot = this->mspScheme->splittableTotal(transId, Splitting::Locations::FIXED);

      // Build a simple balanced split
      ArrayI ids(1);
      ArrayI n0(1);
      ArrayI nN(1);
      SplittingTools::balancedSplit(n0(0), nN(0), tot, this->factor(0), cpuId);

      // Storage for the forward 1D indexes
      std::vector<ArrayI>  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<ArrayI>  bwd1D;
      // Storage for the 2D indexes
      std::vector<ArrayI>  idx2D;
      // Storage for the 3D indexes
      ArrayI  idx3D;

      // Compute the indexes
      ids(0) = cpuId;
      this->mspScheme->fillIndexes(transId, fwd1D, bwd1D, idx2D, idx3D, ids, this->factors(), n0, nN, Splitting::Locations::FIXED);

      // Create TransformResolution object
      return SharedTransformResolution(new TransformResolution(fwd1D, bwd1D, idx2D, idx3D));
   }

   void FixedSplitting::selectGrouper()
   {
      // SINGLE1D or TRANSFORM grouper setup
      #if defined GEOMHDISCC_TRANSGROUPER_SINGLE1D
         this->mGrouper = Splitting::Groupers::SINGLE1D;
      #else
         this->mGrouper = Splitting::Groupers::EQUATION;
      #endif //defined GEOMHDISCC_TRANSGROUPER_SINGLE1D
   }

   int FixedSplitting::computeScore(SharedResolution spResolution)
   {
      // Initialise the score
      double score = 100;

      // Multiply by communication score
      ArrayI comm;
      score *= this->communicationScore(spResolution, comm);

      // Multiply by load balancing score
      Array balance = this->mspScheme->loadWeights();
      score *= this->balancingScore(spResolution, balance);

      // Use additional memory related weighting
      score *= this->mspScheme->memoryScore(spResolution);

      // Select best transform grouper algorithm
      this->selectGrouper();

      return static_cast<int>(score);
   }

}
}
