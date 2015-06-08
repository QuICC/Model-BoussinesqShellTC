/** 
 * @file ParityTools.hpp
 * @brief Implementation of the tools for parity splitting
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PARITYTOOLS_HPP
#define PARITYTOOLS_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Resolutions/Resolution.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   /**
    * @brief Implementation of the tools for parity splitting
    */
   class ParityTools
   {
      public:
         /**
          * @brief Build map of indexes for a regular spatial schemes
          */
         static void splitParity(SharedResolution spRes, const Dimensions::Transform::Id traId, ArrayI& howmany, MatrixI& evenBlocks, MatrixI& oddBlocks);

   };
}
}

#endif // PARITYTOOLS_HPP
