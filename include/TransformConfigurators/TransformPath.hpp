/** 
 * @file TransformPath.hpp
 * @brief This class defines a general transform path
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMPATH_HPP
#define TRANSFORMPATH_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/Arithmetics.hpp"
#include "Enums/FieldIds.hpp"
#include "TransformConfigurators/TransformPathEdge.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class describes a backward transform tree branch
    */
   class TransformPath
   {  
      public:
         /**
          * @brief Contructor for path
          */
         TransformPath(int startId, FieldType::Id fieldId);

         /**
          * @brief Destructor
          */
         ~TransformPath();
         
         /**
          * @brief Add an edge to the tranform path with single ID for output field
          */
         void addEdge(const int projId, const int outId = -1, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Add an edge to the transform path with ID pair for output field
          */
         void addEdge(const int projId, const std::pair<int,int>& outId, Arithmetics::Id arithId);

         /**
          * @brief Get edge
          */
         const TransformPathEdge& edge(const int i) const;

         /**
          * @brief Get the number of edges
          */
         int nEdges() const;

         /**
          * @brief Get starting component ID
          */
         int startId() const;

         /**
          * @brief Get field type ID
          */
         FieldType::Id fieldId() const;

      private:
         /**
          * @brief Starting component for transform path
          */
         int mStartId;

         /**
          * @brief Field type required for transform branch
          */
         FieldType::Id mFieldId;

         /**
          * @brief Vector of edges of the branch
          */
         std::vector<TransformPathEdge> mEdges;
   };

}
}

#endif // TRANSFORMPATH_HPP
