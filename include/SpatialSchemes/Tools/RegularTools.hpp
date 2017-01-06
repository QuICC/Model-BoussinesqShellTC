/** 
 * @file RegularTools.hpp
 * @brief Implementation of the tools for the regular regular schemes
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef REGULARTOOLS_HPP
#define REGULARTOOLS_HPP

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

namespace QuICC {

namespace Schemes {

   /**
    * @brief Implementation of the tools for the regular spatial schemes
    */
   class RegularTools
   {
      public:
         /**
          * @brief Build map of indexes for a regular spatial schemes
          */
         static void buildMap(std::multimap<int,int>& modes, const int i0, const int iN, const ArrayI& j0, const ArrayI& jN, const int c0, const int cN);

         /**
          * @brief Fill the indexes for 2D and 3D for regular spatial schemes
          */
         static void fillIndexes2D3D(std::vector<ArrayI>& idx2D, ArrayI& idx3D, const std::multimap<int,int>& modes);

         /**
          * @brief Fill the indexes for 1D for regular spatial schemes
          */
         static void fillIndexes1D(std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, const ArrayI& idx3D, const int nF1D, const int nB1D);

   };
}
}

#endif // REGULARTOOLS_HPP
