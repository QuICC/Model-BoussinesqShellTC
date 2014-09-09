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

namespace GeoMHDiSCC {

namespace Equations {

   namespace internal
   {
      void addExplicitWrapper(Matrix& rEqField, const int eqStart, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhs);

      template <typename TData> void addExplicitWrapper(TData& rEqField, const int eqStart, const SparseMatrix& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void addExplicitWrapper(MatrixZ& rEqField, const int eqStart, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void addExplicitWrapper(DecoupledZMatrix& rEqField, const int eqStart, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void applyOperatorWrapper(Matrix& rField, const int start, const int rows, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhs);

      void applyOperatorWrapper(MatrixZ& rField, const int start, const int rows, const SparseMatrix& mat, const Eigen::Ref<const MatrixZ>& rhs);

      void applyOperatorWrapper(DecoupledZMatrix& rField, const int start, const int rows, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag);

      void setTopBlock(MatrixZ& rField, const int start, const int rows, const MatrixZ& rhs);

      void setTopBlock(Matrix& rField, const int start, const int rows, const Matrix& rhs);

      void setTopBlock(DecoupledZMatrix& rField, const int start, const int rows, const DecoupledZMatrix& rhs);
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
         virtual void init();

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
          * @brief Compute the source term
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual Datatypes::SpectralScalarType::PointType sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds) = 0;
         
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

         /**
          * \brief Implementation of the coupling definition to python scripts
          */
         void dispatchCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasQI, const bool hasSource, const SharedResolution spRes, const bool allowExplicit = true);

         /**
          * @brief Implementation of base arguments common to all dispatcher
          */
         PyObject* dispatchBaseArguments(const int tupleSize, const ModelOperatorBoundary::Id bcType, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const;

         /**
          * @brief Implementation of model operator dispatcher to python scripts
          */
         void dispatchModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const ModelOperatorBoundary::Id bcType, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const;

         /**
          * @brief Implementation of the galerkin stencil dispatch to python scripts
          */
         void dispatchGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const;

         /**
          * @brief Implementation of the quasi inverse matrix operator dispatch to python scripts
          */
         void dispatchQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const;

         /**
          * @brief Implementation of the explicit linear matrix operator dispatch to python scripts
          */
         void dispatchExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int matIdx, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const;

      private:
         /**
          * @brief Set the galerkin stencil
          */
         virtual void setGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const; // = 0;

         /**
          * @brief Set the quasi inverse matrix operator
          */
         virtual void setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const; // = 0;

         /**
          * @brief Set the explicit linear matrix operator
          */
         virtual void setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int matIdx) const; // = 0;

   };

   /// Typedef for a smart IEquation
   typedef SharedPtrMacro<IEquation> SharedIEquation;

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param eq         Equation
    * @param compId     Equation field component ID
    * @param eqField    Equation field values
    * @param eqStart    Start index for the equation field
    * @param fieldId    Physical field ID
    * @param explicitField   Explicit linear field values
    * @param matIdx     System index
    */
   template <typename TData> void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx);
   template <typename TOperator,typename TData> void computeExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx);

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

      inline void addExplicitWrapper(MatrixZ& rEqField, const int eqStart, const SparseMatrixZ& mat, const Eigen::Ref<const MatrixZ>& rhs)
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

      inline void applyOperatorWrapper(Matrix& rField, const int start, const int rows, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhs)
      {
         int cols = rField.cols();
         rField.block(start, 0, rows, cols) = mat*rhs;
      }

      inline void applyOperatorWrapper(MatrixZ& rField, const int start, const int rows, const SparseMatrix& mat, const Eigen::Ref<const MatrixZ>& rhs)
      {
         int cols = rField.cols();
         rField.block(start, 0, rows, cols) = mat*rhs;
      }

      inline void applyOperatorWrapper(DecoupledZMatrix& rField, const int start, const int rows, const SparseMatrix& mat, const Eigen::Ref<const Matrix>& rhsReal, const Eigen::Ref<const Matrix>& rhsImag)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) = mat*rhsReal;
         rField.imag().block(start, 0, rows, cols) = mat*rhsImag;
      }

      inline void setTopBlock(Matrix& rField, const int start, const int rows, const Matrix& rhs)
      {
         int cols = rField.cols();
         rField.block(start, 0, rows, cols) = rhs.topRows(rows);
      }

      inline void setTopBlock(MatrixZ& rField, const int start, const int rows, const MatrixZ& rhs)
      {
         int cols = rField.cols();
         rField.block(start, 0, rows, cols) = rhs.topRows(rows);
      }

      inline void setTopBlock(DecoupledZMatrix& rField, const int start, const int rows, const DecoupledZMatrix& rhs)
      {
         assert(rField.real().rows() == rField.imag().rows());
         assert(rField.real().cols() == rField.imag().cols());

         int cols = rField.real().cols();
         rField.real().block(start, 0, rows, cols) = rhs.real().topRows(rows);
         rField.imag().block(start, 0, rows, cols) = rhs.imag().topRows(rows);
      }
   }

   template <typename TData> inline void applyQuasiInverse(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& rField, const int start, const int matIdx, const int rhsStart, const TData& rhs)
   {
      // Create pointer to sparse operator
      const SparseMatrix * op = &eq.quasiInverse(compId, matIdx);

      // Get number of rows and cols
      int rows = eq.couplingInfo(compId).galerkinN(matIdx);
      int cols = rField.cols();
      int rhsRows = op->cols();

      internal::applyOperatorWrapper(rField, start, rows, *op, rhs.block(rhsStart, 0, rhsRows, cols));
   }

   template <> inline void applyQuasiInverse<DecoupledZMatrix>(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& rField, const int start, const int matIdx, const int rhsStart, const DecoupledZMatrix& rhs)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      // Create pointer to sparse operator
      const SparseMatrix * op = &eq.quasiInverse(compId, matIdx);

      // Get number of rows and cols
      int rows = eq.couplingInfo(compId).galerkinN(matIdx);
      int cols = rField.real().cols();
      int rhsRows = op->cols();

      internal::applyOperatorWrapper(rField, start, rows, *op, rhs.real().block(rhsStart, 0, rhsRows, cols), rhs.imag().block(rhsStart, 0, rhsRows, cols));
   }

   template <typename TData> inline void applyGalerkinStencil(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& rField, const int start, const int matIdx, const TData& rhs)
   {
      // Create pointer to sparse operator
      const SparseMatrix * op = &eq.galerkinStencil(compId, matIdx);

      // Get number of rows and cols
      int rows = eq.couplingInfo(compId).tauN(matIdx);

      internal::applyOperatorWrapper(rField, 0, rows, *op, rhs.block(start, 0, op->cols(), rhs.cols()));
   }

   template <> inline void applyGalerkinStencil<DecoupledZMatrix>(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& rField, const int start, const int matIdx, const DecoupledZMatrix& rhs)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      // Create pointer to sparse operator
      const SparseMatrix * op = &eq.galerkinStencil(compId, matIdx);

      // Get number of rows and cols
      int rows = eq.couplingInfo(compId).tauN(matIdx);

      internal::applyOperatorWrapper(rField, 0, rows, *op, rhs.real().block(start, 0, op->cols(), rhs.real().cols()), rhs.imag().block(start, 0, op->cols(), rhs.imag().cols()));
   }

   template <typename TData> void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx)
   {
      // Compute with complex linear operator
      if(eq.hasExplicitZLinear(compId, fieldId))
      {
         computeExplicitLinear<SparseMatrixZ>(eq, compId, eqField,  eqStart, fieldId, explicitField, matIdx);
      }

      // Compute with real linear operator
      if(eq.hasExplicitDLinear(compId, fieldId))
      {
         computeExplicitLinear<SparseMatrix>(eq, compId, eqField,  eqStart, fieldId, explicitField, matIdx);
      }
   }

   template <typename TOperator,typename TData> void computeExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx)
   {
      // Create pointer to sparse operator
      const TOperator * op = &eq.explicitLinear<TOperator>(compId, fieldId, matIdx);

      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         /// \mhdBug very bad and slow implementation!
         Eigen::Matrix<Datatypes::SpectralScalarType::PointType,Eigen::Dynamic,1>  tmp(op->cols());
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
         assert(op->cols() == explicitField.slice(mode(0)).rows());

         // Apply operator to field
         internal::addExplicitWrapper(eqField, eqStart, *op, explicitField.slice(mode(0)).col(mode(1)));

      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SINGLE)
      {
         assert(matIdx == 0);

         /// \mhdBug very bad and slow implementation!
         Eigen::Matrix<Datatypes::SpectralScalarType::PointType,Eigen::Dynamic,1>  tmp(op->cols());
         int l = 0;
         for(int k = 0; k < eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); k++)
         {
            for(int j = 0; j < explicitField.slice(k).cols(); j++)
            {
               for(int i = 0; i < explicitField.slice(k).cols(); i++)
               {
                  // Copy slice into flat array
                  tmp(l) = explicitField.point(i,j,k);

                  // increase storage counter
                  l++;
               }
            }
         }

         // Apply operator to field
         internal::addExplicitWrapper(eqField, eqStart, *op, tmp);
      }
   }

   //
   // Dummy specialization
   //

   template <> inline void computeExplicitLinear<SparseMatrixZ,Matrix>(const IEquation& eq, FieldComponents::Spectral::Id compId, Matrix& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx) {};


   
}
}

#endif // IEQUATION_HPP
