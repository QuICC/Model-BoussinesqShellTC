/**
 * @file CouplingInformation.hpp
 * @brief Implemenation of the coupling information of an equation 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef COUPLINGINFORMATION_HPP
#define COUPLINGINFORMATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "Enums/ModelOperator.hpp"
#include "Equations/Tools/IEigenTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implemenation of the coupling information of the equation
    */
   class CouplingInformation
   {
      public:
         /// Typedef to simplify notation for the field data
         typedef std::vector<SpectralFieldId> FieldIdVector;

         /// Typedef for an iterator for the field data
         typedef FieldIdVector::const_iterator  FieldId_iterator;

         /// Typedef for a range iterator for the field coupling data
         typedef std::pair<FieldId_iterator,FieldId_iterator>  FieldId_range;

         /**
          * @brief Enum to specify the type of indexes
          */
         enum IndexType {
            /// Matrix index is slowest index of field
            SLOWEST_SINGLE_RHS = 0,
            /// Matrix index is slowest index of field but has multiple RHS
            SLOWEST_MULTI_RHS,
            /// Matrix index is a mode index
            MODE,
            /// Single matrix (ex. TTT scheme)
            SINGLE,
         };

         /**
          * @brief Enum for the equation type
          */
         enum EquationTypeId {
            /// Equation needs time marching
            PROGNOSTIC = 1,
            /// Equation needs a solver
            DIAGNOSTIC = 2,
            /// Equation is trivial (but with NL, QI, source, explicit)
            TRIVIAL = 3,
            /// Equation does nothing (ie "wrapper" for variable)
            WRAPPER = 4
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
          * @brief Get equation type
          */
         EquationTypeId equationType() const;

         /**
          * @brief Has a nonlinear term?
          */
         bool hasNonlinear() const;

         /**
          * @brief Has a quasi-inverse operator for nonlinear terms?
          */
         bool hasQuasiInverse() const;

         /**
          * @brief Has a source term?
          */
         bool hasSource() const;

         /**
          * @brief Is the system complex?
          */
         bool isComplex() const;

         /**
          * @brief Is Galerkin system?
          */
         bool isGalerkin() const;

         /**
          * @brief Get index type
          */
         IndexType indexType() const;

         /**
          * @brief Get eigen tools
          */
         const IEigenTools& eigenTools() const;

         /**
          * @brief Number of blocks in row of the system
          */
         int nBlocks() const;

         /**
          * @brief Number of systems
          */
         int nSystems() const;

         /**
          * @brief Size of Tau block
          *
          * @param idx System index
          */
         int tauN(const int idx) const;

         /**
          * @brief Size of Galerkin block
          *
          * @param idx System index
          */
         int galerkinN(const int idx) const;

         /**
          * @brief Index shift for the Galerkin block
          *
          * @param dim dimension
          */
         int galerkinShift(const int dim) const;

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
          * @brief Add field to list of implicit timestep fields (coupled solve)
          *
          * @param fieldId Physical ID of the field
          * @param compId  Physical ID of the component
          */
         void addImplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId);

         /**
          * @brief Add field to list of explicit timestep fields
          *
          * @param fieldId Physical ID of the field
          * @param compId  Physical ID of the component
          */
         void addExplicitField(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId, const ModelOperator::Id opId);

         /**
          * @brief Sort the list of implicit fields and set field index accordingly
          *
          * @param fieldId Physical ID of the field
          * @param compId  Physical ID of the component
          */
         void sortImplicitFields(const PhysicalNames::Id fieldId, const FieldComponents::Spectral::Id compId);

         /**
          * @brief Set settings for all systems
          *
          * @param typeId     Type of the solver
          * @param isComplex  Complex flag of solver
          * @param fieldStart Start index of the filed
          */
         void setGeneral(const CouplingInformation::EquationTypeId typeId, const bool isComplex, const int fieldStart);

         /**
          * @brief Set nonlinear flags
          *
          * @param hasNonlinear     Equation requires nonlinear computation?
          * @param hasQuasiInverse  Equation requires quasi-inverse on nonlinear terms?
          */
         void setNonlinear(const bool hasNonlinear, const bool hasQuasiInverse);

         /**
          * @brief Set source flag
          *
          * @param hasSource  Equation requires source term computation?
          */
         void setSource(const bool hasSource);

         /**
          * @brief Set system sizes
          *
          * @param nSystems         Number of systems
          * @param tauNs            Tau block size for each system
          * @param galerkinNs       Galerkin block size for each system
          * @param galerkinShifts   Galerkin index shifts
          * @param rhsCols          Number of columns in RHS for each system
          * @param systemNs         Number of columns in RHS for each system
          */
         void setSizes(const int nSystems, const ArrayI& tauNs, const ArrayI& galerkinNs, const ArrayI& galerkinShifts, const ArrayI& rhsCols, const ArrayI& systemNs);

         /**
          * @brief Set the index type
          *
          * @param id   Dimension type id of index
          */
         void setIndexType(const IndexType id);

         /**
          * @brief set the solver index
          */
         void setSolverIndex(const int idx);

         /**
          * @brief Get iterator to implicit fields
          */
         FieldId_range implicitRange() const;

         /**
          * @brief Get iterator to explicit fields
          */
         FieldId_range explicitRange(const ModelOperator::Id opId) const;

      protected:

      private:
         /**
          * @brief Storage for the implicit fields information
          */
         std::vector<SpectralFieldId>   mImplicitFields;

         /**
          * @brief Storage for the explicit linear fields information
          */
         std::vector<SpectralFieldId>   mExplicitLFields;

         /**
          * @brief Storage for the explicit nonlinear fields information
          */
         std::vector<SpectralFieldId>   mExplicitNLFields;

         /**
          * @brief Storage for the explicit nextstep fields information
          */
         std::vector<SpectralFieldId>   mExplicitNSFields;

         /**
          * @brief Storage for the equation type
          */
         EquationTypeId mEquationType;

         /**
          * @brief Storage for the nonlinear flag
          */
         bool mHasNonlinear;

         /**
          * @brief Storage for the quasi-inverse flag
          */
         bool mHasQuasiInverse;

         /**
          * @brief Storage for the source flag
          */
         bool mHasSource;

         /**
          * @brief Storage for the complex flag
          */
         bool mIsComplex;

         /**
          * @brief Storage for the galerkin flag
          */
         bool mIsGalerkin;

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
          * @brief Shared IEigenTools
          */
         SharedIEigenTools mspEigenTools;

         /**
          * @brief Storage for the Tau block sizes
          */
         ArrayI mTauNs;

         /**
          * @brief Storage for the Galerkin block sizes
          */
         ArrayI mGalerkinNs;

         /**
          * @brief Storage for the Galerkin index shifts
          */
         ArrayI mGalerkinShifts;

         /**
          * @brief Storage for the number of RHS in solver
          */
         ArrayI mRhsCols;

         /**
          * @brief Storage for the system sizes
          */
         ArrayI mSystemNs;
   };
}
}

#endif // COUPLINGINFORMATION_HPP
