/**
 * @file IEquation.hpp
 * @brief Base building block for the implementation of an equation 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IEQUATION_HPP
#define IEQUATION_HPP

// First include
//
#include "Python/PythonHeader.hpp"

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
#include "Base/MatrixOperationsInternal.hpp"
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "Enums/ModelOperator.hpp"
#include "Enums/ModelOperatorBoundary.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/CouplingInformation.hpp"
#include "Equations/EquationData.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Variables/VariableRequirement.hpp"
#include "Simulation/SimulationBoundary.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base building block for the implementation of an equation
    */
   class IEquation : public EquationData
   {
      public:
         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IEquation();

         /**
          * @brief Access the shared resolution
          */
         virtual SharedResolution spRes() const = 0;

         /**
          * @brief Initialise the equation
          */
         virtual void init(const SharedSimulationBoundary spBcIds);

         /**
          * @brief Generic model operator dispatcher to python scripts
          */
         virtual void buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const ModelOperatorBoundary::Id bcType) const; // = 0;

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Use the nonlinear computation as physical field values
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void useNonlinear(const Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id);

         /**
          * @brief Compute the source term
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual Datatypes::SpectralScalarType::PointType sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Compute the boundary value
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual Datatypes::SpectralScalarType::PointType boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Initialise the spectral equation matrices
          */
         virtual void initSpectralMatrices() = 0;

         /**
          * @brief Implementation of the galerkin stencil dispatch to python scripts
          */
         void dispatchGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx, const SharedResolution spRes, const std::vector<MHDFloat>& eigs, const bool makeSquare = false) const;
         
      protected:
         /**
          * @brief Set the equation variable requirements
          */
         virtual void setRequirements() = 0;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() = 0;

         /**
          * @brief Set the default nonlinear components
          */
         virtual void setNLComponents() = 0;

         /**
          * @brief Initialise the spectral equation matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initSpectralMatricesComponent(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

         /**
          * \brief Implementation of the coupling definition to python scripts
          */
         void dispatchCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasSource, const SharedResolution spRes, const bool hasBoundaryValue = false, const bool allowExplicit = true);

         /**
          * @brief Implementation of base arguments common to all dispatcher
          */
         PyObject* dispatchBaseArguments(const int tupleSize, const ModelOperatorBoundary::Id bcType, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const;

         /**
          * @brief Implementation of model operator dispatcher to python scripts
          */
         void dispatchModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const ModelOperatorBoundary::Id bcType, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const;

         /**
          * @brief Implementation of the explicit matrix operator dispatch to python scripts
          */
         void dispatchExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const ModelOperator::Id opId, const SpectralFieldId fieldId, const int matIdx, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const;

         /**
          * @brief Update the stored value with the solver solution (real data)
          */
         virtual MHDFloat updateStoredSolution(const MHDFloat newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k);

         /**
          * @brief Update the stored value with the solver solution (complex data)
          */
         virtual MHDComplex updateStoredSolution(const MHDComplex newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k);

         /**
          * @brief Tune solution obtained from solver
          */
         virtual void tuneSolution(const FieldComponents::Spectral::Id compId, const int k);

      private:
         /**
          * @brief Initialise the quasi-inverse spectral  matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initQIMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

         /**
          * @brief Initialise the explicit spectral  matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          * @param opId    Type of explicit operator
          */
         void initExplicitMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId, const ModelOperator::Id opId);

         /**
          * @brief Set the galerkin stencil
          */
         virtual void setGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const; // = 0;

         /**
          * @brief Set the explicit matrix operator
          */
         virtual void setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const ModelOperator::Id opId, const SpectralFieldId fieldId, const int matIdx) const; // = 0;

   };

   /// Typedef for a smart IEquation
   typedef SharedPtrMacro<IEquation> SharedIEquation;

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param eq            Equation
    * @param compId        Equation field component ID
    * @param rSolverField  Solver field values
    * @param eqStart       Start index for the equation field
    * @param fieldId       Physical field ID
    * @param explicitField Explicit linear field values
    * @param matIdx        System index
    */
   template <typename TData> void addExplicitTerm(const IEquation& eq, const ModelOperator::Id opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx);
   template <typename TOperator,typename TData> void computeExplicitTerm(const IEquation& eq, const ModelOperator::Id opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx);

   /**
    * @brief Apply the quasi-inverse operator
    *
    * @param eq         Equation
    * @param compId     Equation field component ID
    * @param rField     Output field
    * @param start      Start index in linear storage
    * @param matIdx     System index
    * @param rhsStart   Start index in RHS data
    * @param rhs        RHS field data
    */
   template <typename TData> void applyQuasiInverse(const IEquation& eq, TData& rField, const int start, const int matIdx, const int rhsStart, const TData& rhs);
   template <> void applyQuasiInverse<DecoupledZMatrix>(const IEquation& eq, DecoupledZMatrix& rField, const int start, const int matIdx, const int rhsStart, const DecoupledZMatrix& rhs);

   /**
    * @brief Apply the galerkin stencil operator
    *
    * @param eq         Equation
    * @param compId     Equation field component ID
    * @param rField     Output field
    * @param start      Start index in linear storage
    * @param matIdx     System index
    * @param rhs        RHS field data
    */
   template <typename TData> void applyGalerkinStencil(const IEquation& eq, TData& rField, const int start, const int matIdx, const TData& rhs);
   template <> void applyGalerkinStencil<DecoupledZMatrix>(const IEquation& eq, DecoupledZMatrix& rField, const int start, const int matIdx, const DecoupledZMatrix& rhs);

   inline MHDFloat IEquation::updateStoredSolution(const MHDFloat newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k)
   {
      return newData;
   }

   inline MHDComplex IEquation::updateStoredSolution(const MHDComplex newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k)
   {
      return newData;
   }

   template <typename TData> inline void applyQuasiInverse(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& rField, const int start, const int matIdx, const int rhsStart, const TData& rhs)
   {
      if(eq.hasQID(compId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.quasiInverse<SparseMatrix>(compId, matIdx);

         // Get number of rows and cols
         int cols = rField.cols();
         int rhsRows = op->cols();

         Datatypes::internal::addMatrixProduct(rField, start, *op, rhs.block(rhsStart, 0, rhsRows, cols));

      } else if(eq.hasQIZ(compId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.quasiInverse<SparseMatrixZ>(compId, matIdx);

         // Get number of rows and cols
         int cols = rField.cols();
         int rhsRows = op->cols();

         Datatypes::internal::addMatrixProduct(rField, start, *op, rhs.block(rhsStart, 0, rhsRows, cols));
      }
   }

   template <> inline void applyQuasiInverse<DecoupledZMatrix>(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& rField, const int start, const int matIdx, const int rhsStart, const DecoupledZMatrix& rhs)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      if(eq.hasQID(compId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.quasiInverse<SparseMatrix>(compId, matIdx);

         // Get number of rows and cols
         int cols = rField.real().cols();
         int rhsRows = op->cols();

         Datatypes::internal::addMatrixProduct(rField, start, *op, rhs.real().block(rhsStart, 0, rhsRows, cols), rhs.imag().block(rhsStart, 0, rhsRows, cols));

      } else if(eq.hasQIZ(compId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.quasiInverse<SparseMatrixZ>(compId, matIdx);

         // Get number of rows and cols
         int cols = rField.real().cols();
         int rhsRows = op->cols();

         Datatypes::internal::addMatrixProduct(rField, start, *op, rhs.real().block(rhsStart, 0, rhsRows, cols), rhs.imag().block(rhsStart, 0, rhsRows, cols));
      }
   }

   template <typename TData> inline void applyGalerkinStencil(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& rField, const int start, const int matIdx, const TData& rhs)
   {
      // Create pointer to sparse operator
      const SparseMatrix * op = &eq.galerkinStencil(compId, matIdx);

      Datatypes::internal::setMatrixProduct(rField, 0, *op, rhs.block(start, 0, op->cols(), rhs.cols()));
   }

   template <> inline void applyGalerkinStencil<DecoupledZMatrix>(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& rField, const int start, const int matIdx, const DecoupledZMatrix& rhs)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      // Create pointer to sparse operator
      const SparseMatrix * op = &eq.galerkinStencil(compId, matIdx);

      Datatypes::internal::setMatrixProduct(rField, 0, *op, rhs.real().block(start, 0, op->cols(), rhs.real().cols()), rhs.imag().block(start, 0, op->cols(), rhs.imag().cols()));
   }

   template <typename TData> void addExplicitTerm(const IEquation& eq, const ModelOperator::Id opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx)
   {
      // Compute with complex linear operator
      if(eq.hasExplicitZTerm(opId, compId, fieldId))
      {
         computeExplicitTerm<SparseMatrixZ>(eq, opId, compId, rSolverField,  eqStart, fieldId, explicitField, matIdx);
      }

      // Compute with real linear operator
      if(eq.hasExplicitDTerm(opId, compId, fieldId))
      {
         computeExplicitTerm<SparseMatrix>(eq, opId, compId, rSolverField,  eqStart, fieldId, explicitField, matIdx);
      }
   }

   template <typename TOperator,typename TData> void computeExplicitTerm(const IEquation& eq, const ModelOperator::Id opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx)
   {
      // Create pointer to sparse operator
      const TOperator * op = &eq.explicitOperator<TOperator>(opId, compId, fieldId, matIdx);

      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
      {
         Eigen::Matrix<Datatypes::SpectralScalarType::PointType,Eigen::Dynamic,1>  tmp(op->cols());
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            // Initialise storage to zero
            tmp.setZero();
            int l;
            int j_;
            int dimI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
               int corrDim = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
            #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
            for(int j = 0; j < explicitField.slice(matIdx).cols(); j++)
            {
               j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
               #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                  j_ -= corrDim;
               #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
               for(int i = 0; i < explicitField.slice(matIdx).rows(); i++)
               {
                  // Compute correct position
                  l = j_ + i;

                  // Copy field value into storage
                  tmp(l) = explicitField.point(i,j,matIdx);
               }
            }
         #else
            int k = 0;
            for(int j = 0; j < explicitField.slice(matIdx).cols(); j++)
            {
               for(int i = 0; i < explicitField.slice(matIdx).rows(); i++)
               {
                  // Copy slice into flat array
                  tmp(k) = explicitField.point(i,j,matIdx);

                  // increase storage counter
                  k++;
               }
            }
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

         // Apply operator to field
         Datatypes::internal::addMatrixProduct(rSolverField, eqStart, *op, tmp);

      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
      {
         // Apply operator to field
         Datatypes::internal::addMatrixProduct(rSolverField, eqStart, *op, explicitField.slice(matIdx));

      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         // Get mode indexes
         ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

         // Assert correct sizes
         assert(op->cols() == explicitField.slice(mode(0)).rows());

         // Apply operator to field
         Datatypes::internal::addMatrixProduct(rSolverField, eqStart, *op, explicitField.slice(mode(0)).col(mode(1)));

      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SINGLE)
      {
         assert(matIdx == 0);

         /// \mhdBug very bad and slow implementation!
         Eigen::Matrix<Datatypes::SpectralScalarType::PointType,Eigen::Dynamic,1>  tmp(op->cols());
         int l = 0;
         int k_;
         int j_;
         #ifdef QUICC_SPATIALDIMENSION_3D
            int dimK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
         #else
            int dimK = 1;
         #endif //QUICC_SPATIALDIMENSION_3D
         int dimJ = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         for(int k = 0; k < eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); k++)
         {
            k_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)*dimK;
            for(int j = 0; j < explicitField.slice(k).cols(); j++)
            {
               j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
               for(int i = 0; i < explicitField.slice(k).rows(); i++)
               {
                  // Compute correct position
                  l = k_ + j_ + i;

                  // Copy slice into flat array
                  tmp(l) = explicitField.point(i,j,k);
               }
            }
         }

         // Apply operator to field
         Datatypes::internal::addMatrixProduct(rSolverField, eqStart, *op, tmp);
      }
   }

   //
   // Dummy specialization
   //

   template <> inline void computeExplicitTerm<SparseMatrixZ,Matrix>(const IEquation& eq, const ModelOperator::Id opId, FieldComponents::Spectral::Id compId, Matrix& rSolverField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx) {};


   
}
}

#endif // IEQUATION_HPP
