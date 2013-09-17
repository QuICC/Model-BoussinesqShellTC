/** 
 * @file SplittingDescription.cpp
 * @brief Source of the base of the implementation of the load splitting algorithm description
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "LoadSplitter/Algorithms/SplittingDescription.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {

   SplittingDescription::SplittingDescription(const Splitting::Algorithms::Id algorithm, const Splitting::Groupers::Id grouper, const int dims, const ArrayI& factors, const MHDFloat score)
      :  algorithm(algorithm), grouper(grouper), dims(dims), factors(factors), score(score)
   {
   }

   SplittingDescription::~SplittingDescription()
   {
   }

}
}
