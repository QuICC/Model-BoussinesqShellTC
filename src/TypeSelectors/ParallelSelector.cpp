/** 
 * @file ParallelSelector.cpp
 * @brief Definitions of the ParallelSelector grouper selection
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

   #ifdef QUICC_MPI
      #ifdef QUICC_TRANSGROUPER_SINGLE1D
      } else if(descr.grouper == Splitting::Groupers::SINGLE1D)
      {
         setGrouper<Splitting::Groupers::SINGLE1D>(descr.algorithm, spFwdGrouper, spBwdGrouper);
      #endif //QUICC_TRANSGROUPER_SINGLE1D
      #ifdef QUICC_TRANSGROUPER_SINGLE2D
      } else if(descr.grouper == Splitting::Groupers::SINGLE2D)
      {
         setGrouper<Splitting::Groupers::SINGLE2D>(descr.algorithm, spFwdGrouper, spBwdGrouper);
      #endif //QUICC_TRANSGROUPER_SINGLE2D
      #ifdef QUICC_TRANSGROUPER_TRANSFORM
      } else if(descr.grouper == Splitting::Groupers::TRANSFORM)
      {
         setGrouper<Splitting::Groupers::TRANSFORM>(descr.algorithm, spFwdGrouper, spBwdGrouper);
      #endif //QUICC_TRANSGROUPER_TRANSFORM
   #endif //QUICC_MPI
      } else
      {
         throw Exception("Unknown transform grouper setup");
      }
   }
}
}
