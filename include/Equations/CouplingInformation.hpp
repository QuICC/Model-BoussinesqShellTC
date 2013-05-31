/** \file CouplingInformation.hpp
 *  \brief Implemenation of the coupling information of an equation
 */

#ifndef COUPLINGINFORMATION_HPP
#define COUPLINGINFORMATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implemenation of the coupling information of the equation
    */
   class CouplingInformation
   {
      public:
         /// Typedef to simplify notation for the field data
         typedef std::vector<SpectralFieldId> FieldIdVector;

         /// Typedef for an iterator for the field data
         typedef FieldIdVector::const_iterator  field_iterator;

         /// Typedef for a range iterator for the field coupling data
         typedef std::pair<field_iterator,field_iterator>  field_iterator_range;

         /**
          * @brief Enum to specify the type of indexes
          */
         enum IndexType {
            /// Matrix index is slowest index of field
            SLOWEST = 0,
            /// Matrix index is a mode index
            MODE
         };

         /**
          * @brief Simple constructor
          */
         CouplingInformation();

         /**
          * @brief Simple empty destructor
          */
         virtual ~CouplingInformation();

         /**
          * @brief Is the system complex?
          */
         bool isComplex() const;

         /**
          * @brief Get index type
          */
         IndexType indexType() const;

         /**
          * @brief Number of blocks in row of the system
          */
         int nBlocks() const;

         /**
          * @brief Number of systems
          */
         int nSystems() const;

         /**
          * @brief Size of block
          *
          * @param idx System index
          */
         int blockN(const int idx) const;

         /**
          * @brief Size of system
          *
          * @param idx System index
          */
         int systemN(const int idx) const;

         /**
          * @brief Index of the field within the system
          */
         int fieldIndex() const;

         /**
          * @brief Index of solver
          */
         int solverIndex() const;

         /**
          * @brief Start index in field values
          */
         int fieldStart() const;

         /**
          * @brief Number of RHS columns in solve
          *
          * @param idx System index
          */
         int rhsCols(const int idx) const;

         /**
          * @brief Number of explicit linear fields
          */
         int nExplicit() const;

         /**
          * @brief Add field to list of implicit timestep fields (coupled solve)
          *
          * @param fieldId Physical ID of the field
          * @param compId  Physical ID of the component
          * @param isSelf  Is field equation field?
          */
         void addImplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId, const bool isSelf);

         /**
          * @brief Add field to list of explicit timestep fields
          *
          * @param fieldId Physical ID of the field
          * @param compId  Physical ID of the component
          */
         void addExplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId);

         /**
          * @brief Set settings for all systems
          *
          * @param solverIndex   Index of the solver
          * @param isComplex     Complex flag of solver
          * @param fieldStart    Start index of the filed
          */
         void setGeneral(const int solverIndex, const bool isComplex, const int fieldStart);

         /**
          * @brief Set system sizes
          *
          * @param nSystems   Number of systems
          * @param blockNs    Block size for each system
          * @param rhsCols    Number of columns in RHS for each system
          */
         void setSizes(const int nSystems, const ArrayI& blockNs, const ArrayI& rhsCols);

         /**
          * @brief Set the index type
          *
          * @param id   Dimension type id of index
          */
         void setIndexType(const IndexType id);

         /**
          * @brief Get iterator to implicit fields
          */
         field_iterator_range implicitRange() const;

         /**
          * @brief Get iterator to explicit fields
          */
         field_iterator_range explicitRange() const;

      protected:

      private:
         /**
          * @brief Storage for the implicit fields information
          */
         std::vector<SpectralFieldId>   mImplicitFields;

         /**
          * @brief Storage for the explicit fields information
          */
         std::vector<SpectralFieldId>   mExplicitFields;

         /**
          * @brief Storage for the complex flag
          */
         bool mIsComplex;

         /**
          * @brief Storage for the index type
          */
         IndexType mIndexType;

         /**
          * @brief Storage for the number of systems
          */
         int mNSystems;

         /**
          * @brief Storage for the field index
          */
         int mFieldIndex;

         /**
          * @brief Storage for the solver index
          */
         int mSolverIndex;

         /**
          * @brief Storage for the field start index
          */
         int mFieldStart;

         /**
          * @brief Storage for the block sizes
          */
         ArrayI mBlockNs;

         /**
          * @brief Storage for the number of RHS in solver
          */
         ArrayI mRhsCols;
   };
}
}

#endif // COUPLINGINFORMATION_HPP
