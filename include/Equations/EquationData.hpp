/**
 * @file EquationData.hpp
 * @brief Lowest building block for the implementation of an equation 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONDATA_HPP
#define EQUATIONDATA_HPP

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
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/CouplingInformation.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Variables/VariableRequirement.hpp"
#include "Simulation/SimulationBoundary.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Lowest building block for the implementation of an equation
    */
   class EquationData
   {
      public:
         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit EquationData(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~EquationData();

         /**
          * @brief Get name ID of the unknown
          */
         PhysicalNames::Id name() const;

         /**
          * @brief Get scalar variable
          *
          * @param name Physical name of the field
          */
         const Datatypes::ScalarVariableType& scalar(PhysicalNames::Id name) const;

         /**
          * @brief Get vector variable
          *
          * @param name Physical name of the field
          */
         const Datatypes::VectorVariableType& vector(PhysicalNames::Id name) const;

         /**
          * @brief Set the smart pointer to the scalar field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the scalar field
          */
         void setField(PhysicalNames::Id name, Datatypes::SharedScalarVariableType spField);

         /**
          * @brief Set the smart pointer to the vector field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         void setField(PhysicalNames::Id name, Datatypes::SharedVectorVariableType spField);

         /**
          * @brief Get the galerkin stencil matrix
          *
          * @param compId  Field component ID
          * @param j       Matrix index
          */
         const SparseMatrix& galerkinStencil(const FieldComponents::Spectral::Id compId, const int j) const;

         /**
          * @brief Get the quasi-inverse matrix for the nonlinear terms
          *
          * @param compId  Field component ID
          * @param j       Matrix index
          */
         const SparseMatrix& quasiInverse(const FieldComponents::Spectral::Id compId, const int j) const;

         /**
          * @brief Check if real explicit linear matrices exist
          *
          * @param compId  Field component ID
          * @param fieldId Spectral field ID
          */
         bool hasExplicitDLinear(const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const;

         /**
          * @brief Check if complex explicit linear matrices exist
          *
          * @param compId  Field component ID
          * @param fieldId Spectral field ID
          */
         bool hasExplicitZLinear(const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const;

         /**
          * @brief Get the explicit linear matrices
          *
          * @param compId  Field component ID
          * @param fieldId Spectral field ID
          * @param j       Matrix index
          */
         template <typename TOperator> const TOperator& explicitLinear(const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const int j) const;

         /**
          * @brief Get the coupling information
          *
          * @param compId  Field component ID
          */
         const CouplingInformation&  couplingInfo(const FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Get map of field storage requirements information
          *
          * \mhdBug Ultimatively this should depend on component
          */
         const VariableRequirement& requirements() const;

         /**
          * @brief Get map of field storage requirements information
          */
         const FieldRequirement& requirements(PhysicalNames::Id id) const;

         /**
          * @brief Get the equation parameters
          */
         const EquationParameters& eqParams() const;

         /**
          * @brief Get the list of boundary conditions
          */
         const SimulationBoundary& bcIds() const;

         /**
          * @brief Set the solver index for the coupling information
          */
         void setSolverIndex(const FieldComponents::Spectral::Id, const int idx);

      protected:
         /**
          * @brief Set the unknown name of equation
          */
         void setName(PhysicalNames::Id name);

         /**
          * @brief Set scalar variable
          *
          * @param name Physical name of the field
          */
         Datatypes::ScalarVariableType& rScalar(PhysicalNames::Id name);

         /**
          * @brief Set vector variable
          *
          * @param name Physical name of the field
          */
         Datatypes::VectorVariableType& rVector(PhysicalNames::Id name);

         /**
          * @brief Storage for the variable requirements
          *
          * \mhdBug Ultimatively this should depend on component
          */
         VariableRequirement mRequirements;

         /**
          * @brief Coupling information of the equation
          */
         std::map<FieldComponents::Spectral::Id, CouplingInformation>  mCouplingInfos;

         /**
          * @brief Map of component and explicit linear matrices (real operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> > mLDMatrices;

         /**
          * @brief Map of component and explicit linear matrices (complex operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> > mLZMatrices;

         /**
          * @brief Map of component and nonlinear term multiplication matrices
          */
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> > mNLMatrices;

         /**
          * @brief Map of component and galerkin stencil matrices
          */
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> > mGStencils;

         /**
          * @brief Storage for the shared boundary condition list
          */
         SharedSimulationBoundary mspBcIds;

      private:
         /**
          * @brief Storage for smart equation parameters
          */
         SharedEquationParameters   mspEqParams;

         /**
          * @brief Map of name and pointer for the scalar variables
          */
         std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>  mScalars;

         /**
          * @brief Map of name and pointer for the vector variables
          */
         std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>  mVectors;

         /**
          * @brief Name ID of the unknown
          */
         PhysicalNames::Id mName;
   };

   /// Typedef for a smart EquationData
   typedef SharedPtrMacro<EquationData> SharedEquationData;
}
}

#endif // EQUATIONDATA_HPP
