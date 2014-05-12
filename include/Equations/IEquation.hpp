/**
 * @file IEquation.hpp
 * @brief Base building block for the implementation of an equation 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IEQUATION_HPP
#define IEQUATION_HPP

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
#include "Enums/ModelOperator.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/CouplingInformation.hpp"
#include "Equations/EquationData.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/BoundaryMethodSelector.hpp"
#include "Variables/VariableRequirement.hpp"
#include "Simulation/SimulationBoundary.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   namespace internal
   {
      void addExplicitWrapper(Matrix& rEqField, const int eqStart, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhs);

      template <typename TData> void addExplicitWrapper(TData& rEqField, const int eqStart, const SparseMatrix& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void addExplicitWrapper(MatrixZ& rEqField, const int eqStart, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void addExplicitWrapper(DecoupledZMatrix& rEqField, const int eqStart, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs);

      template <typename TData> void applyQuasiInverseWrapper(TData& rField, const int start, const int rows, const SparseMatrix& mat);

      template <> void applyQuasiInverseWrapper<DecoupledZMatrix>(DecoupledZMatrix& rField, const int start, const int rows, const SparseMatrix& mat);
   }

   /**
    * @brief Base building block for the implementation of an equation
    */
   class IEquation : public EquationData
   {
      public:
         /**
          * @brief Simple constructor
          *
          * \param pyName     Python script name
          * \param spEqParams Shared equation parameters
          */
         explicit IEquation(const std::string& pyName, SharedEquationParameters spEqParams);

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
         virtual void init();

         /**
          * @brief Generic model operator dispatcher to python scripts
          */
         virtual void buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const bool hasBoundary) const; // = 0;

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Compute the source term
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual MHDComplex sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds) = 0;

         /**
          * @brief Get the boundary condition coordinator
          *
          * @param compId  Component ID
          */
         virtual const Boundary::CoordinatorSelector& bcCoord(FieldComponents::Spectral::Id compId) const = 0;

         /**
          * @brief Set the boundary condition coordinator
          *
          * @param compId  Component ID
          */
         virtual Boundary::CoordinatorSelector& rBcCoord(FieldComponents::Spectral::Id compId) = 0;
         
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
          * @brief Initialise the spectral equation matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initSpectralMatricesComponent(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

      private:
         /**
          * @brief Set the quasi inverse matrix operator
          */
         void setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix &mat) const;

         /**
          * @brief Set the explicit linear matrix operator
          */
         void setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const;

   };

   /// Typedef for a smart IEquation
   typedef SharedPtrMacro<IEquation> SharedIEquation;

   /**
    * @brief Apply quasi-inverse to the nonlinear term
    *
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void applyQuasiInverse(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param compId     Equation field component ID
    * @param eqField    Equation field values
    * @param eqStart    Start index for the equation field
    * @param fieldId    Physical field ID
    * @param explicitField   Explicit linear field values
    * @param matIdx     System index
    */
   template <typename TData> void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx);

   namespace internal
   {
      inline void addExplicitWrapper(Matrix& rEqField, const int eqStart, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhs)
      {
         rEqField.block(eqStart, 0, mat.rows(), rEqField.cols()) += mat*rhs;
      }

      template <typename TData> inline void addExplicitWrapper(TData& rEqField, const int eqStart, const SparseMatrix& mat, const Eigen::Ref<const MatrixZ>& rhs)
      {
         int rows = mat.rows();
         int cols = rEqField.real().cols();
         rEqField.real().block(eqStart, 0, rows, cols) += mat*rhs.real();
         rEqField.imag().block(eqStart, 0, rows, cols) += mat*rhs.imag();
      }

      inline void addExplicitWrapper(MatrixZ& rEqField, const int eqStart, const SparseMatrixZ& mat, const typename Eigen::Ref<const MatrixZ>& rhs)
      {
         rEqField.block(eqStart, 0, mat.rows(), rEqField.cols()) += mat*rhs;
      }

      inline void addExplicitWrapper(DecoupledZMatrix& rEqField, const int eqStart, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs)
      {
         assert(rEqField.real().rows() == rEqField.imag().rows());
         assert(rEqField.real().cols() == rEqField.imag().cols());

         int rows = mat.rows();
         int cols = rEqField.real().cols();
         rEqField.real().block(eqStart, 0, rows, cols) += mat.real()*rhs.real() - mat.imag()*rhs.imag();
         rEqField.imag().block(eqStart, 0, rows, cols) += mat.real()*rhs.imag() + mat.imag()*rhs.real();
      }

      template <typename TData> inline void applyQuasiInverseWrapper(TData& rField, const int start, const int rows, const SparseMatrix& mat)
      {
         int cols = rField.cols();
         rField.block(start, 0, rows, cols) = mat*rField.block(start, 0, rows, cols);
      }

      template <> inline void applyQuasiInverseWrapper<DecoupledZMatrix>(DecoupledZMatrix& rField, const int start, const int rows, const SparseMatrix& mat)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) = mat*rField.real().block(start, 0, rows, cols);
         rField.imag().block(start, 0, rows, cols) = mat*rField.imag().block(start, 0, rows, cols);
      }
   }

   template <typename TData> void applyQuasiInverse(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Apply quasi inverse
      if(eq.couplingInfo(compId).hasQuasiInverse())
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.quasiInverse(compId, matIdx);

         // Get number of rows
         int rows = eq.couplingInfo(compId).blockN(matIdx);

         // Safety asserts
         assert(op->rows() == op->cols());
         assert(op->cols() == rows);

         // Multiply nonlinear term by quasi-inverse
         internal::applyQuasiInverseWrapper(storage, start, rows, *op);
      }
   }

   template <typename TData> void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx)
   {
      // Compute with complex linear operator
      if(eq.hasExplicitZLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.explicitZLinear(compId, fieldId, matIdx);

         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            /// \mhdBug very bad and slow implementation!
            Eigen::Matrix<Datatypes::SpectralScalarType::PointType,Eigen::Dynamic,1>  tmp(op->rows());
            int k = 0;
            for(int j = 0; j < explicitField.slice(matIdx).cols(); j++)
            {
               for(int i = 0; i < explicitField.slice(matIdx).cols(); i++)
               {
                  // Copy slice into flat array
                  tmp(k) = explicitField.point(i,j,matIdx);

                  // increase storage counter
                  k++;
               }
            }

            // Apply operator to field
            internal::addExplicitWrapper(eqField, eqStart, *op, tmp);

         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            // Get mode indexes
            ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

            // Assert correct sizes
            assert(op->rows() == explicitField.slice(mode(0)).rows());

            // Apply operator to field
            internal::addExplicitWrapper(eqField, eqStart, *op, explicitField.slice(mode(0)).col(mode(1)));
         }
      }

      // Compute with real linear operator
      if(eq.hasExplicitDLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.explicitDLinear(compId, fieldId, matIdx);

         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            /// \mhdBug very bad and slow implementation!
            Eigen::Matrix<Datatypes::SpectralScalarType::PointType,Eigen::Dynamic,1>  tmp(op->rows());
            int k = 0;
            for(int j = 0; j < explicitField.slice(matIdx).cols(); j++)
            {
               for(int i = 0; i < explicitField.slice(matIdx).cols(); i++)
               {
                  // Copy slice into flat array
                  tmp(k) = explicitField.point(i,j,matIdx);

                  // increase storage counter
                  k++;
               }
            }

            // Apply operator to field
            internal::addExplicitWrapper(eqField, eqStart, *op, tmp);

         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            // Get mode indexes
            ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

            // Assert correct sizes
            assert(op->rows() == explicitField.slice(mode(0)).rows());

            // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
            internal::addExplicitWrapper(eqField, eqStart, *op, (*op)*explicitField.slice(mode(0)).col(mode(1)));
         }
      }
   }
   
}
}

#endif // IEQUATION_HPP
