/** 
 * @file TransformTree.cpp
 * @brief Source of the implementation of the tranform tree
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
#include "TransformConfigurators/TransformPathEdge.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   TransformPathEdge::TransformPathEdge(const int opId, const int outId, const Arithmetics::Id arithId)
      :mOpId(opId), mArithId(arithId)
   {
      this->mOutId.push_back(outId);
   }

   TransformPathEdge::TransformPathEdge(const int opId, const std::pair<int,int>& outId, const Arithmetics::Id arithId)
      :mOpId(opId), mArithId(arithId)
   {
      this->mOutId.push_back(outId.first);
      this->mOutId.push_back(outId.second);
   }

   TransformPathEdge::~TransformPathEdge()
   {
   }

   int TransformPathEdge::opId() const
   {
      return this->mOpId;
   }

   const std::vector<int>& TransformPathEdge::outId() const
   {
      return this->mOutId;
   }

   Arithmetics::Id TransformPathEdge::arithId() const
   {
      return this->mArithId;
   }



}
}
