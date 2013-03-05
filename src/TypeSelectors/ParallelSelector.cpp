/** \file ParallelSelector.cpp
 *  \brief Definitions of the ParallelSelector grouper selection
 */

// System includes
//

// External includes
//

// Class include
//
#include "TypeSelectors/ParallelSelector.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {
   
   void setGrouper(const SplittingDescription& descr, Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper)
   {
      if(descr.grouper == Splitting::Groupers::EQUATION)
      {
         setGrouper<Splitting::Groupers::EQUATION>(descr.algorithm, spFwdGrouper, spBwdGrouper);
      } else if(descr.grouper == Splitting::Groupers::SINGLE1D)
      {
         setGrouper<Splitting::Groupers::SINGLE1D>(descr.algorithm, spFwdGrouper, spBwdGrouper);
      } else if(descr.grouper == Splitting::Groupers::SINGLE2D)
      {
         setGrouper<Splitting::Groupers::SINGLE2D>(descr.algorithm, spFwdGrouper, spBwdGrouper);
      } else if(descr.grouper == Splitting::Groupers::TRANSFORM)
      {
         setGrouper<Splitting::Groupers::TRANSFORM>(descr.algorithm, spFwdGrouper, spBwdGrouper);
      }
   }
}
}
