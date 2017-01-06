/** 
 * @file ParityTools.cpp
 * @brief Source of the tools for parity splitting
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <set>
#include <map>
#include <tr1/tuple>

// External includes
//

// Class include
//
#include "SpatialSchemes/Tools/ParityTools.hpp"

// Project includes
//

namespace QuICC {

namespace Schemes {

   void ParityTools::splitParityL(SharedResolution spRes, const Dimensions::Transform::Id traId, ArrayI& howmany, MatrixI& evenBlocks, MatrixI& oddBlocks)
   {
      // Get number of transforms
      std::vector<std::tr1::tuple<int,int,int> >::iterator it;
      std::vector<std::tr1::tuple<int,int,int> > even;
      std::vector<std::tr1::tuple<int,int,int> > odd;
      int idx = 0;
      int previous = -1;
      for(int i = 0; i < spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         if(spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT3D>(i)%2 == 0)
         {
            if(previous == 0)
            {
               std::tr1::get<1>(even.back()) += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);
            } else
            {
               even.push_back(std::tr1::make_tuple(idx, spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i), spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT3D>(i)));
            }

            idx += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);

            previous = 0;
         } else
         {
            if(previous == 1)
            {
               std::tr1::get<1>(odd.back()) += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);
            } else
            {
               odd.push_back(std::tr1::make_tuple(idx, spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i),spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT3D>(i)));
            }

            idx += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);

            previous = 1;
         }
      }

      evenBlocks.resize(even.size(), 3);
      oddBlocks.resize(odd.size(), 3);

      int i = 0;
      for(it = even.begin(); it != even.end(); ++it, ++i)
      {
         evenBlocks(i,0) = std::tr1::get<0>(*it);
         evenBlocks(i,1) = std::tr1::get<1>(*it);
         evenBlocks(i,2) = std::tr1::get<2>(*it);
      }

      i = 0;
      for(it = odd.begin(); it != odd.end(); ++it, ++i)
      {
         oddBlocks(i,0) = std::tr1::get<0>(*it);
         oddBlocks(i,1) = std::tr1::get<1>(*it);
         oddBlocks(i,2) = std::tr1::get<2>(*it);
      }

      howmany.resize(2);
      howmany(0) = evenBlocks.col(1).sum();
      howmany(1) = oddBlocks.col(1).sum();
   }

   void ParityTools::splitParityM(SharedResolution spRes, const Dimensions::Transform::Id traId, ArrayI& howmany, MatrixI& evenBlocks, MatrixI& oddBlocks)
   {
      // Get number of transforms
      std::vector<std::tr1::tuple<int,int,int> >::iterator it;
      std::vector<std::tr1::tuple<int,int,int> > even;
      std::vector<std::tr1::tuple<int,int,int> > odd;
      int idx = 0;
      int previous = -1;
      for(int i = 0; i < spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         for(int j = 0; j < spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i); j++)
         {
            if(spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT2D>(j,i)%2 == 0)
            {
               if(previous == 0)
               {
                  std::tr1::get<1>(even.back()) += 1;
               } else
               {
                  even.push_back(std::tr1::make_tuple(idx, 1, spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT2D>(j,i)));
               }

               idx += 1;

               previous = 0;
            }
            else
            {
               if(previous == 1)
               {
                  std::tr1::get<1>(odd.back()) += 1;
               } else
               {
                  odd.push_back(std::tr1::make_tuple(idx, 1, spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT2D>(j,i)));
               }

               idx += 1;

               previous = 1;
            }
         }
      }

      evenBlocks.resize(even.size(), 3);
      oddBlocks.resize(odd.size(), 3);

      int i = 0;
      for(it = even.begin(); it != even.end(); ++it, ++i)
      {
         evenBlocks(i,0) = std::tr1::get<0>(*it);
         evenBlocks(i,1) = std::tr1::get<1>(*it);
         evenBlocks(i,2) = std::tr1::get<2>(*it);
      }

      i = 0;
      for(it = odd.begin(); it != odd.end(); ++it, ++i)
      {
         oddBlocks(i,0) = std::tr1::get<0>(*it);
         oddBlocks(i,1) = std::tr1::get<1>(*it);
         oddBlocks(i,2) = std::tr1::get<2>(*it);
      }

      howmany.resize(2);
      howmany(0) = evenBlocks.col(1).sum();
      howmany(1) = oddBlocks.col(1).sum();
   }

}
}
