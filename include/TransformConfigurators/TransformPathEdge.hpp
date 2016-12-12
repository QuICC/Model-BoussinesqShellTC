/** 
 * @file TransformPathEdge.hpp
 * @brief This class defines an edge in transform path
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMPATHEDGE_HPP
#define TRANSFORMPATHEDGE_HPP

// Configuration includes
// 

// System includes
//
#include <set>
#include <vector>

// External includes
//

// Project includes
//
#include "Enums/Arithmetics.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines an edge in a transform paht
    */
   class TransformPathEdge
   {  
      public:
         /**
          * @brief Contructor for edge with single ID for output field
          */
         TransformPathEdge(const int opId, const int outId, const Arithmetics::Id arithId);

         /**
          * @brief Contructor for edge with ID pair for output field
          */
         TransformPathEdge(const int opId, const std::pair<int,int>& outId, const Arithmetics::Id arithId);

         /**
          * @brief Destructor
          */
         ~TransformPathEdge();

         /**
          * @brief Get ID of operator
          */
         int opId() const;

         /**
          * @brief Get IDs of output field
          */
         const std::vector<int>& outId() const;

         /**
          * @brief Get ID of arithmetic operation
          */
         Arithmetics::Id arithId() const;

      private:
         /**
          * @brief ID of operatr
          */
         int mOpId;

         /**
          * @brief ID of the arithmetic operation
          */
         Arithmetics::Id   mArithId;

         /**
          * @brief IDs of the output field
          */
         std::vector<int>  mOutId;
   };

}
}

#endif // TRANSFORMPATHEDGE_HPP
