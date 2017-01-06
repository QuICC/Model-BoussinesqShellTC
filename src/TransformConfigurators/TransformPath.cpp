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
#include "TransformConfigurators/TransformPath.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

   TransformPath::TransformPath(int startId, FieldType::Id fieldId)
      :mStartId(startId), mFieldId(fieldId)
   {
   }

   TransformPath::~TransformPath()
   {
   }
         
   void TransformPath::addEdge(const int opId, const int outId, Arithmetics::Id arithId)
   {
      this->mEdges.push_back(TransformPathEdge(opId, outId, arithId));
   }

   void TransformPath::addEdge(const int opId, const std::pair<int,int>& outId, Arithmetics::Id arithId)
   {
      this->mEdges.push_back(TransformPathEdge(opId, outId, arithId));
   }

   const TransformPathEdge& TransformPath::edge(const int i) const
   {
      return this->mEdges.at(i);
   }

   int TransformPath::nEdges() const
   {
      return this->mEdges.size();
   }

   int TransformPath::startId() const
   {
      return this->mStartId;
   }

   FieldType::Id TransformPath::fieldId() const
   {
      return this->mFieldId;
   }

}
}
