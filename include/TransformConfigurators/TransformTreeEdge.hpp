/** 
 * @file TransformTreeEdge.hpp
 * @brief This class defines a general transform tree edge
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMTREEEDGE_HPP
#define TRANSFORMTREEEDGE_HPP

// Configuration includes
// 

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Enums/Arithmetics.hpp"
#include "Enums/FieldIds.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   // Forward declaration
   class TransformTreeEdge;

   /**
    * @brief Specialisation to end recursion
    */
   class TransformTreeEdge
   {
      public:
         /// Useful typedefs
         typedef std::vector<TransformTreeEdge>::const_iterator EdgeType_iterator;
         typedef std::pair<EdgeType_iterator,EdgeType_iterator> EdgeType_range;

         /**
          * @brief Contructor for operation
          */
         TransformTreeEdge(const int op, const int n);

         /**
          * @brief Destructor
          */
         ~TransformTreeEdge();

         /**
          * @brief Get operator ID
          */
         template <typename TId> TId opId() const;

         /**
          * @brief Get the number of edges
          */
         int nEdges(const int depth) const;

         /**
          * @brief Get vector ranges of edges
          */
         EdgeType_range edgeRange() const;

         /**
          * @brief Add a edge to tree
          */
         TransformTreeEdge& addEdge(const int op, const int n);

         /**
          * @brief Get end component
          */
          template <typename TId> TId endId(const int i = 0) const;

          /**
           * @brief Get the field type
           */
          FieldType::Id fieldId() const;

          /**
           * @brief Get the arithmetic operation
           */
          Arithmetics::Id arithId() const;

         /**
          * @brief Set end component as vector of IDs
          */
         void setEnd(const std::vector<int>& id, const FieldType::Id type, const Arithmetics::Id arith);

         /**
          * @brief Set end component as a single ID
          */
         void setEnd(const int id, const FieldType::Id type, const Arithmetics::Id arith);

      private:
         /**
          * @brief Operation attached to edge
          */
         int  mOpId;

         /**
          * @brief Multiplicity of edge (used by how many branches)
          */
         int mN;

         /**
          * @brief Field type
          */
         FieldType::Id  mFieldId;

         /**
          * @brief Arithmetic operation to store result
          */
         Arithmetics::Id   mArithId;

         /**
          * @brief Physical output component
          */
         std::vector<int> mEndId;

         /**
          * @brief Vector of connected egdes
          */
         std::vector<TransformTreeEdge>  mEdges;
   };

   template <typename TId> inline TId TransformTreeEdge::opId() const
   {
      return static_cast<TId>(this->mOpId);
   }

   template <typename TId> inline TId TransformTreeEdge::endId(const int i) const
   {
      return static_cast<TId>(this->mEndId.at(i));
   }
}
}

#endif // TRANSFORMTREEEDGE_HPP
