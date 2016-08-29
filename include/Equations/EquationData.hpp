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
#include <limits>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "Enums/SolveTiming.hpp"
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
          * @brief Check if real quasi inverse matrices exist
          *
          * @param compId  Field component ID
          */
         bool hasQID(const FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Check if complex quasi inverse matrices exist
          *
          * @param compId  Field component ID
          */
         bool hasQIZ(const FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Check if real explicit matrices exist
          *
          * @param compId  Field component ID
          * @param fieldId Spectral field ID
          */
         bool hasExplicitDTerm(const ModelOperator::Id opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const;

         /**
          * @brief Check if complex explicit matrices exist
          *
          * @param compId  Field component ID
          * @param fieldId Spectral field ID
          */
         bool hasExplicitZTerm(const ModelOperator::Id opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const;

         /**
          * @brief Get the quasi inverse matrices
          *
          * @param compId  Field component ID
          * @param j       Matrix index
          */
         template <typename TOperator> const TOperator& quasiInverse(const FieldComponents::Spectral::Id compId, const int j) const;

         /**
          * @brief Get the explicit matrices
          *
          * @param compId  Field component ID
          * @param fieldId Spectral field ID
          * @param j       Matrix index
          */
         template <typename TOperator> const TOperator& explicitOperator(const ModelOperator::Id opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const int j) const;

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
          * @brief Get map of imposed field storage requirements information
          *
          * \mhdBug Ultimatively this should depend on component
          */
         const VariableRequirement& imposedRequirements() const;

         /**
          * @brief Get map of field storage requirements information
          */
         const FieldRequirement& requirements(PhysicalNames::Id id) const;

         /**
          * @brief Get map of imposed field storage requirements information
          */
         const FieldRequirement& imposedRequirements(PhysicalNames::Id id) const;

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

         /**
          * @brief Timing of the solver for the equation
          */
         SolveTiming::Id  solveTiming() const;

         /**
          * @brief Current simulation time to allow for timedependent implementations
          */
         MHDFloat  time() const;

         /**
          * @brief Set current simulation time to allow for timedependent implementations
          *
          * @param time       Current simulation time
          * @param finished   Flag for completed multistage timestep
          */
         virtual void  setTime(const MHDFloat time, const bool finished);

         /**
          * @brief Get the nonlinear integration components order
          */
         const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& nlComponents() const;

      protected:
         /**
          * @brief Set the unknown name of equation
          */
         void setName(PhysicalNames::Id name);

         /**
          * @brief Set solver timing
          */
         void setSolveTiming(const SolveTiming::Id time);

         /**
          * @brief Add a nonlinear integration component
          */
         void addNLComponent(const FieldComponents::Spectral::Id compId, const int flag);

         /**
          * @brief Update field requirements information
          */
         FieldRequirement& updateFieldRequirements(PhysicalNames::Id id);

         /**
          * @brief Set map of component and explicit matrices (real operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> >& rEDMatrices(const ModelOperator::Id opId);

         /**
          * @brief Set map of component and explicit matrices (complex operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> >& rEZMatrices(const ModelOperator::Id opId);

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
          * @brief Storage for the imposed variable requirements
          *
          * \mhdBug Ultimatively this should depend on component
          */
         VariableRequirement mImposedRequirements;

         /**
          * @brief Coupling information of the equation
          */
         std::map<FieldComponents::Spectral::Id, CouplingInformation>  mCouplingInfos;

         /**
          * @brief Map of component and quasi inverse matrices (real operators)
          */
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> > mQIDMatrices;

         /**
          * @brief Map of component and quasi inverse matrices (complex operators)
          */
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrixZ> > mQIZMatrices;

         /**
          * @brief Map of component and galerkin stencil matrices
          */
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> > mGStencils;

         /**
          * @brief Storage for the shared boundary condition list
          */
         SharedSimulationBoundary mspBcIds;

         /**
          * @brief Storage for the solve timing
          */
         SolveTiming::Id   mSolveTiming;

         /**
          * @brief Nonlinear integration component order
          */
         std::vector<std::pair<FieldComponents::Spectral::Id,int> >   mNLComponents;

      private:
         /**
          * @brief Storage for smart equation parameters
          */
         SharedEquationParameters   mspEqParams;

         /**
          * @brief Name ID of the unknown
          */
         MHDFloat mTime;

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

         /**
          * @brief Map of component and explicit linear matrices (real operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> > mELDMatrices;

         /**
          * @brief Map of component and explicit linear matrices (complex operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> > mELZMatrices;

         /**
          * @brief Map of component and explicit nonlinear matrices (real operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> > mENLDMatrices;

         /**
          * @brief Map of component and explicit nonlinear matrices (complex operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> > mENLZMatrices;

         /**
          * @brief Map of component and explicit nextstep matrices (real operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> > mENSDMatrices;

         /**
          * @brief Map of component and explicit nextstep matrices (complex operators)
          */
         std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> > mENSZMatrices;
   };

   /// Typedef for a smart EquationData
   typedef SharedPtrMacro<EquationData> SharedEquationData;

   /**
    * @brief Update time average
    */
   template <typename TData> TData incrementTimeAverage(const TData avg, const TData newData, const MHDFloat time, const MHDFloat timestep);
   MHDFloat incrementTimeAverage(const MHDComplex avg, const MHDFloat newData, const MHDFloat time, const MHDFloat timestep);

   /**
    * @brief Don't update time average
    */
   template <typename TData> TData noupdateTimeAverage(const TData avg, const TData newData);
   MHDFloat noupdateTimeAverage(const MHDComplex avg, const MHDFloat newData);


   template <typename TData> TData incrementTimeAverage(const TData avg, const TData newData, const MHDFloat time, const MHDFloat timestep)
   {
      MHDFloat stepWeight = timestep/time;
      MHDFloat avgWeight = (time-timestep)/time;

      return avgWeight*avg + stepWeight*newData;
   }

   template <typename TData> TData noupdateTimeAverage(const TData avg, const TData newData)
   {
      return avg;
   }
}
}

#endif // EQUATIONDATA_HPP
