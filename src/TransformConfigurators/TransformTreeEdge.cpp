/** 
 * @file TransformTreeEdge.cpp
 * @brief Source of the implementation of the tranform tree edge
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
#include "TransformConfigurators/TransformTree.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   TransformTreeEdge::TransformTreeEdge(const int op, const int n)
      :mOpId(op), mN(n), mFieldId(FieldType::SCALAR), mArithId(Arithmetics::SET)
   {
   }

   TransformTreeEdge::~TransformTreeEdge()
   {
   }

   int TransformTreeEdge::nEdges(const int depth) const
   {
      if(depth == 0)
      {
         return this->mEdges.size();

      } else
      {
         int n = 0;
         for(std::vector<TransformTreeEdge>::const_iterator it = this->mEdges.begin(); it != this->mEdges.end(); ++it)
         {
            n += it->nEdges(depth - 1);
         }

         return n;
      }
   }

   TransformTreeEdge::EdgeType_range TransformTreeEdge::edgeRange() const
   {
      return std::make_pair(this->mEdges.begin(), this->mEdges.end());
   }

   TransformTreeEdge& TransformTreeEdge::addEdge(const int op, const int n)
   {
      this->mEdges.push_back(TransformTreeEdge(op, n));

      return this->mEdges.back();
   }

   FieldType::Id TransformTreeEdge::fieldId() const
   {
      return this->mFieldId;
   }

   Arithmetics::Id TransformTreeEdge::arithId() const
   {
      return this->mArithId;
   }

   void TransformTreeEdge::setEnd(const int id, const FieldType::Id type, const Arithmetics::Id arith)
   {
      this->mEndId.push_back(id);

      this->mFieldId = type;

      this->mArithId = arith;
   }

   void TransformTreeEdge::setEnd(const std::vector<int>& id, const FieldType::Id type, const Arithmetics::Id arith)
   {
      this->mEndId = id;

      this->mFieldId = type;

      this->mArithId = arith;
   }

}
}
