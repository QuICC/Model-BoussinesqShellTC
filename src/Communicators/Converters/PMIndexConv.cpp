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
#include "Enums/DimensionTools.hpp"

namespace QuICC {

namespace Parallel {

   PMIndexConv::PMIndexConv(SharedResolution spRes, const Dimensions::Transform::Id id)
      : mspTResFwd(spRes->cpu()->dim(id)), mspTResBwd(spRes->cpu()->dim(Dimensions::jump(id,1))), mSimN(spRes->sim()->dim(static_cast<Dimensions::Simulation::Id>(static_cast<int>(id)+1), Dimensions::Space::SPECTRAL))
   {
      /// \mhdBug Calculation of mSimN with double static_cast is not ideal
   }

   PMIndexConv::~PMIndexConv()
   {
   }

}
}
