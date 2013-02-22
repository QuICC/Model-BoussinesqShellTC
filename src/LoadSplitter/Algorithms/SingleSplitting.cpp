/** \file SingleSplitting.cpp
 *  \brief Source of the implementation of a single load splitting algorithm
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
#include "LoadSplitter/Algorithms/SingleSplitting.hpp"

// Project includes
//
#include "Exception/Exception.hpp"
#include "Base/IO/Ascii/FormatToolbox.hpp"

namespace GeoMHDiSCC {

   SingleSplitting::SingleSplitting(const int id, const int nCpu, const ArrayI& dim, Splittings::Locations::Id split)
      : SplittingAlgorithm(id, nCpu, dim, Splittings::Algorithms::SINGLE1D), mSplit(split)
   {
      // Make sure the right algorithm is set
      if(split == Splittings::Locations::SECOND)
      {
         this->mAlgo = Splittings::Algorithms::SINGLE2D;
      }

      // Factorise N_cpu
      this->factoriseNCpu(1);

      // Filter factors
      this->filterFactors();
   }

   bool SingleSplitting::applicable() const
   {
      bool status = true;

      // Only works with at least 2 dimensions
      status = (status && (this->dims() > 1));

      // Make sure split and dimensions are compatible
      if(this->dims() > 2 || this->mSplit == Splittings::Locations::FIRST)
      {
         // All CPUs should have something to do
         int tot;
         for(int i = 0; i < this->dims(); i++)
         {
            tot = this->mspScheme->splittableTotal(i, this->mSplit);

            status = (status && (tot >= this->factor(0)));
         }
      } else
      {
         status = false;
      }

      // Check for scheme specific conditions
      status = (status && this->mspScheme->applicable());

      return status;
   }

   SharedTransformResolution  SingleSplitting::splitDimension(const int dim, const int id)
   {
      // Get size of the splittable dimension(s)
      int tot = this->mspScheme->splittableTotal(dim, this->mSplit);

      // Build a simple balanced split
      ArrayI ids(1);
      ArrayI n0(1);
      ArrayI nN(1);
      this->balancedSplit(n0(0), nN(0), tot, this->factor(0), id);

      // Storage for the forward 1D indexes
      std::vector<ArrayI>  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<ArrayI>  bwd1D;
      // Storage for the 2D indexes
      std::vector<ArrayI>  idx2D;
      // Storage for the 3D indexes
      ArrayI  idx3D;

      // Compute the indexes
      ids(0) = id;
      this->mspScheme->fillIndexes(dim, fwd1D, bwd1D, idx2D, idx3D, ids, this->factors(), n0, nN, this->mSplit);

      // Create TransformResolution object
      return SharedTransformResolution(new TransformResolution(fwd1D, bwd1D, idx2D, idx3D));
   }

   void SingleSplitting::selectGrouper()
   {
      // Different splitting directions require a different treatment
      if(this->mSplit == Splittings::Locations::FIRST)
      {
         // SINGLE1D or TRANSFORM grouper setup
         #if defined GEOMHDISCC_TRANSGROUPER_SINGLE1D || defined GEOMHDISCC_TRANSGROUPER_TRANSFORM
            this->mGrouper = Splittings::Groupers::SINGLE1D;
         #else
            this->mGrouper = Splittings::Groupers::EQUATION;
         #endif //defined GEOMHDISCC_TRANSGROUPER_SINGLE1D || defined GEOMHDISCC_TRANSGROUPER_TRANSFORM
      } else
      {
         // SINGLE2D or TRANSFORM grouper setup
         #if defined GEOMHDISCC_TRANSGROUPER_SINGLE2D || defined GEOMHDISCC_TRANSGROUPER_TRANSFORM
            this->mGrouper = Splittings::Groupers::SINGLE2D;
         #else
            this->mGrouper = Splittings::Groupers::EQUATION;
         #endif //defined GEOMHDISCC_TRANSGROUPER_SINGLE2D || defined GEOMHDISCC_TRANSGROUPER_TRANSFORM
      }
   }

   int SingleSplitting::computeScore(SharedResolution spResolution)
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
