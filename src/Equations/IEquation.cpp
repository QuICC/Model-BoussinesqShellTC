/** \file IEquation.cpp
 *  \brief Source of building block for the implementation of a time dependend evolution equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IEquation::IEquation(SharedEquationParameters spEqParams)
      : EquationData(spEqParams)
   {
   }

   IEquation::~IEquation()
   {
   }

   void IEquation::init()
   {
      this->setCoupling();
   }

   void IEquation::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // This implementation should never get called!
      throw Exception("Activated nonlinear term without implementation!");
   }

   MHDComplex IEquation::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      // This implementation should never get called!
      throw Exception("Activated source term without implementation!");

      return MHDComplex();
   }

   DecoupledZSparse  IEquation::operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id comp, const int matIdx) const
   {
      throw Exception("operatorRow: dummy implementation was called!");
   }

   void applyQuasiInverse(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start)
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
         storage.first.block(start, 0, rows, storage.first.cols()) = (*op)*storage.first.block(start, 0, rows, storage.first.cols());
         storage.second.block(start, 0, rows, storage.second.cols()) = (*op)*storage.second.block(start, 0, rows, storage.second.cols());
      }
   }

   void applyQuasiInverse(const IEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start)
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
         storage.block(start, 0, rows, storage.cols()) = (*op)*storage.block(start, 0, rows, storage.cols());
      }
   }

   void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& eqField, const int eqStart, SpectralFieldId fieldId, const DecoupledZMatrix& linField, const int linStart, const int matIdx)
   {
      // Safety asserts
      assert(eqField.first.cols() == linField.first.cols());
      assert(eqField.second.cols() == linField.second.cols());

      // Compute with complex linear operator
      if(eq.hasExplicitZLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.explicitZLinear(compId, fieldId, matIdx);

         // Get number of rows and columns of block
         int rows = op->rows();
         int cols = eqField.first.cols();

         // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
         eqField.first.block(eqStart, 0, rows, cols) += op->real()*linField.first.block(linStart, 0, rows, cols);
         eqField.first.block(eqStart, 0, rows, cols) -= op->imag()*linField.second.block(linStart, 0, rows, cols);
         eqField.second.block(eqStart, 0, rows, cols) += op->real()*linField.second.block(linStart, 0, rows, cols);
         eqField.second.block(eqStart, 0, rows, cols) += op->imag()*linField.first.block(linStart, 0, rows, cols);
      }

      // Compute with real linear operator
      if(eq.hasExplicitDLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.explicitDLinear(compId, fieldId, matIdx);

         // Get number of rows and columns of block
         int rows = op->rows();
         int cols = eqField.first.cols();

         // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
         eqField.first.block(eqStart, 0, rows, cols) += (*op)*linField.first.block(linStart, 0, rows, cols);
         eqField.second.block(eqStart, 0, rows, cols) += (*op)*linField.second.block(linStart, 0, rows, cols);
      }
   }

   void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& eqField, const int eqStart, SpectralFieldId fieldId, const MatrixZ& linField, const int linStart, const int matIdx)
   {
      // Safety asserts
      assert(eqField.first.cols() == linField.cols());
      assert(eqField.second.cols() == linField.cols());

      if(eq.hasExplicitZLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.explicitZLinear(compId, fieldId, matIdx);

         // Get number of rows and columns of block
         int rows = op->rows();
         int cols = eqField.first.cols();

         // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
         eqField.first.block(eqStart, 0, rows, cols) += op->real()*linField.block(linStart, 0, rows,cols).real();
         eqField.first.block(eqStart, 0, rows, cols) -= op->imag()*linField.block(linStart, 0, rows, cols).imag();
         eqField.second.block(eqStart, 0, rows, cols) += op->real()*linField.block(linStart, 0, rows, cols).imag();
         eqField.second.block(eqStart, 0, rows, cols) += op->imag()*linField.block(linStart, 0, rows, cols).real();
      }

      if(eq.hasExplicitDLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.explicitDLinear(compId, fieldId, matIdx);

         // Get number of rows and columns of block
         int rows = op->rows();
         int cols = eqField.first.cols();

         // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
         eqField.first.block(eqStart, 0, rows, cols) += (*op)*linField.block(linStart, 0, rows, cols).real();
         eqField.second.block(eqStart, 0, rows, cols) += (*op)*linField.block(linStart, 0, rows, cols).imag();
      }
   }

   void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& eqField, const int eqStart, SpectralFieldId fieldId, const DecoupledZMatrix& linField, const int linStart, const int matIdx)
   {
      // Safety asserts
      assert(eqField.cols() == linField.first.cols());
      assert(eqField.cols() == linField.second.cols());

      if(eq.hasExplicitZLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.explicitZLinear(compId, fieldId, matIdx);

         // Get number of rows and columns of block
         int rows = op->rows();
         int cols = eqField.cols();

         // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
         eqField.block(eqStart, 0, rows, cols).real() += op->real()*linField.first.block(linStart, 0, rows, cols);
         eqField.block(eqStart, 0, rows, cols).real() -= op->imag()*linField.second.block(linStart, 0, rows, cols);
         eqField.block(eqStart, 0, rows, cols).imag() += op->imag()*linField.second.block(linStart, 0, rows, cols);
         eqField.block(eqStart, 0, rows, cols).imag() += op->imag()*linField.first.block(linStart, 0, rows, cols);
      }

      if(eq.hasExplicitDLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.explicitDLinear(compId, fieldId, matIdx);

         // Get number of rows and columns of block
         int rows = op->rows();
         int cols = eqField.cols();

         // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
         eqField.block(eqStart, 0, rows, cols).real() += (*op)*linField.first.block(linStart, 0, rows, cols);
         eqField.block(eqStart, 0, rows, cols).imag() += (*op)*linField.second.block(linStart, 0, rows, cols);
      }
   }

   void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& eqField, const int eqStart, SpectralFieldId fieldId, const MatrixZ& linField, const int linStart, const int matIdx)
   {
      // Safety asserts
      assert(eqField.cols() == linField.cols());

      if(eq.hasExplicitZLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.explicitZLinear(compId, fieldId, matIdx);

         // Get number of rows and columns of block
         int rows = op->rows();
         int cols = eqField.cols();

         // Apply operator to field: eq += op*lin
         eqField.block(eqStart, 0, rows, cols) += (*op)*linField.block(linStart, 0, rows, cols);
      }

      if(eq.hasExplicitDLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.explicitDLinear(compId, fieldId, matIdx);

         // Get number of rows and columns of block
         int rows = op->rows();
         int cols = eqField.cols();

         eqField.block(eqStart, 0, rows, cols).real() += (*op)*linField.block(linStart, 0, rows, cols).real();
         eqField.block(eqStart, 0, rows, cols).imag() += (*op)*linField.block(linStart, 0, rows, cols).imag();
      }
   }

   void quasiInverseBlock(const IEquation& eq, SparseMatrix& mat)
   {
      throw Exception("quasiInverseBlock: dummy implementation called!");
   }

   void timeBlock(const IEquation& eq, DecoupledZSparse& mat, const MHDFloat k)
   {
      throw Exception("timeBlock: dummy implementation called!");
   }

   void linearBlock(const IEquation& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      throw Exception("linearBlock: dummy implementation called!");
   }

   void boundaryBlock(const IEquation& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      throw Exception("boundaryBlock: dummy implementation called!");
   }
}
}
