/** 
 * @file TransformEdge.hpp
 * @brief This template describes a transform edge and its connection further down the tree
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMEDGE_HPP
#define TRANSFORMEDGE_HPP

// Configuration includes
// 

// System includes
//
#include <vector>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   // Forward declaration
   template <typename TOpA, typename TOpB, typename TOpC> class TransformEdge;

   /**
    * @brief This template describes a transform edge and its connection further down the tree
    */
   template <typename TOpA, typename TOpB, typename TOpC> class TransformEdge
   {  
      public:
         /// Typedefs to simply notations
         typedef TransformEdge<TOpB,TOpC,void> EdgeType;
         typedef typename std::vector<EdgeType>::const_iterator EdgeType_iterator;
         typedef std::pair<EdgeType_iterator,EdgeType_iterator> EdgeType_range;

         /**
          * @brief Contructor for operation
          */
         TransformEdge(const TOpA op, const int n);

         /**
          * @brief Destructor
          */
         ~TransformEdge();

         /**
          * @brief Get projector operator ID
          */
         TOpA projector() const;

         /**
          * @brief Get the number of edges
          */
         int nEdges() const;

         /**
          * @brief Get vector ranges of edges
          */
         EdgeType_range edgeRange() const;

         /**
          * @brief Add a edge to tree
          */
         EdgeType& addEdge(const TOpB, const int n);

      private:
         /**
          * @brief Operation attached to edge
          */
         TOpA  mOpId;

         /**
          * @brief Multiplicity of edge (used by how many branches)
          */
         int mN;

         /**
          * @brief Vector of connected egdes
          */
         std::vector<EdgeType>  mEdges;
   };

   /**
    * @brief Specialisation to end recursion
    */
   template <typename TOpA> class TransformEdge<TOpA,void,void>
   {
      public:
         /**
          * @brief Contructor for operation
          */
         TransformEdge(const TOpA op, const int n);

         /**
          * @brief Destructor
          */
         ~TransformEdge();

         /**
          * @brief Get projector operator ID
          */
         TOpA projector() const;

         /**
          * @brief Get Physical component
          */
          FieldComponents::Physical::Id physId() const;

          /**
           * @brief Get the field type
           */
          FieldType::Id fieldId() const;

         /**
          * @brief Set Physical component
          */
         void setPhysical(const FieldComponents::Physical::Id id, const FieldType::Id type);

      private:
         /**
          * @brief Operation attached to edge
          */
         TOpA  mOpId;

         /**
          * @brief Multiplicity of edge (used by how many branches)
          */
         int mN;

         /**
          * @brief Field type
          */
         FieldType::Id  mFieldId;

         /**
          * @brief Physical output component
          */
         FieldComponents::Physical::Id mPhysId;
   };


   template <typename TOpA, typename TOpB, typename TOpC> TransformEdge<TOpA,TOpB,TOpC>::TransformEdge(const TOpA op, const int n)
      :mOpId(op), mN(n)
   {
   }

   template <typename TOpA> TransformEdge<TOpA,void,void>::TransformEdge(const TOpA op, const int n)
      :mOpId(op), mN(n)
   {
   }

   template <typename TOpA, typename TOpB, typename TOpC> TransformEdge<TOpA,TOpB,TOpC>::~TransformEdge()
   {
   }

   template <typename TOpA> TransformEdge<TOpA,void,void>::~TransformEdge()
   {
   }

   template <typename TOpA, typename TOpB, typename TOpC> TOpA TransformEdge<TOpA,TOpB,TOpC>::projector() const
   {
      return this->mOpId;
   }

   template <typename TOpA, typename TOpB, typename TOpC> int TransformEdge<TOpA,TOpB,TOpC>::nEdges() const
   {
      return this->mEdges.size();
   }

   template <typename TOpA, typename TOpB, typename TOpC> typename TransformEdge<TOpA,TOpB,TOpC>::EdgeType_range TransformEdge<TOpA,TOpB,TOpC>::edgeRange() const
   {
      return std::make_pair(this->mEdges.begin(), this->mEdges.end());
   }

   template <typename TOpA, typename TOpB, typename TOpC> TransformEdge<TOpB,TOpC,void>& TransformEdge<TOpA,TOpB,TOpC>::addEdge(const TOpB op, const int n)
   {
      this->mEdges.push_back(TransformEdge<TOpB,TOpC,void>(op, n));

      return this->mEdges.back();
   }

   template <typename TOpA> TOpA TransformEdge<TOpA,void,void>::projector() const
   {
      return this->mOpId;
   }

   template <typename TOpA> FieldComponents::Physical::Id TransformEdge<TOpA,void,void>::physId() const
   {
      return this->mPhysId;
   }

   template <typename TOpA> FieldType::Id TransformEdge<TOpA,void,void>::fieldId() const
   {
      return this->mFieldId;
   }

   template <typename TOpA> void TransformEdge<TOpA,void,void>::setPhysical(const FieldComponents::Physical::Id id, const FieldType::Id type)
   {
      this->mPhysId = id;

      this->mFieldId = type;
   }

}
}

#endif // TRANSFORMEDGE_HPP
