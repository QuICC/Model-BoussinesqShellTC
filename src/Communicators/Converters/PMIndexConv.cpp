/** 
 * @file PMIndexConv.cpp
 * @brief Source of the index converter that splits array for positive and negative FFT frequencies
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Communicators/Converters/PMIndexConv.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {

   PMIndexConv::PMIndexConv(SharedCTransformResolution spTRes)
      : mspTRes(spTRes)
   {
   }

   PMIndexConv::~PMIndexConv()
   {
   }

}
}
