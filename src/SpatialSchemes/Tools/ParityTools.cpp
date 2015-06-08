/** 
 * @file ParityTools.cpp
 * @brief Source of the tools for parity splitting
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/Tools/ParityTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Schemes {

   void ParityTools::splitParity(SharedResolution spRes, const Dimensions::Transform::Id traId, ArrayI& howmany, MatrixI& evenBlocks, MatrixI& oddBlocks)
   {
      // Get number of transforms
      std::vector<std::pair<int,int> >::iterator it;
      std::vector<std::pair<int,int> > even;
      std::vector<std::pair<int,int> > odd;
      int idx = 0;
      int previous = -1;
      for(int i = 0; i < spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         if(spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT3D>(i)%2 == 0)
         {
            if(previous == 0)
            {
               even.back().second += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);
            } else
            {
               even.push_back(std::make_pair(idx, spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i)));
            }

            idx += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);

            previous = 0;
         } else
         {
            if(previous == 1)
            {
               odd.back().second += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);
            } else
            {
               odd.push_back(std::make_pair(idx, spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i)));
            }

            idx += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);

            previous = 1;
         }
      }

      evenBlocks.resize(even.size(), 2);
      oddBlocks.resize(odd.size(), 2);

      int i = 0;
      for(it = even.begin(); it != even.end(); ++it, ++i)
      {
         evenBlocks(i,0) = it->first;
         evenBlocks(i,1) = it->second;
      }

      i = 0;
      for(it = odd.begin(); it != odd.end(); ++it, ++i)
      {
         oddBlocks(i,0) = it->first;
         oddBlocks(i,1) = it->second;
      }

      howmany.resize(2);
      howmany(0) = evenBlocks.col(1).sum();
      howmany(1) = oddBlocks.col(1).sum();
   }

}
}
