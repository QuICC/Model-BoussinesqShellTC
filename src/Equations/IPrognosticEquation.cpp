/** \file IPrognosticEquation.cpp
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
#include "Equations/IPrognosticEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IPrognosticEquation::IPrognosticEquation(SharedEquationParameters spEqParams)
      : EquationData(spEqParams)
   {
   }

   IPrognosticEquation::~IPrognosticEquation()
   {
   }

   void IPrognosticEquation::init()
   {
      this->setCoupling();
   }

   void IPrognosticEquation::computeLinear(FieldComponents::Spectral::Id compId, DecoupledZMatrix& eqField, const int eqStart, SpectralFieldId fieldId, const DecoupledZMatrix& linField, const int linStart, const int matIdx) const
   {
      // Safety asserts
      assert(eqField.first.cols() == linField.first.cols());
      assert(eqField.second.cols() == linField.second.cols());

      // Create key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(compId, fieldId);

      if(this->mLZMatrices.count(key) > 0)
      {
         // Get number of rows and columns of block
         int rows = this->mLZMatrices.find(key)->second.at(matIdx).rows();
         int cols = eqField.first.cols();

         // Create pointer to sparse operator
         const SparseMatrixZ * op = &this->mLZMatrices.find(key)->second.at(matIdx);

         // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
         eqField.first.block(eqStart, 0, rows, cols) += op->real()*linField.first.block(linStart, 0, rows, cols);
         eqField.first.block(eqStart, 0, rows, cols) -= op->imag()*linField.second.block(linStart, 0, rows, cols);
         eqField.second.block(eqStart, 0, rows, cols) += op->real()*linField.second.block(linStart, 0, rows, cols);
         eqField.second.block(eqStart, 0, rows, cols) += op->imag()*linField.first.block(linStart, 0, rows, cols);
      }

      if(this->mLDMatrices.count(key) > 0)
      {
         // Get number of rows and columns of block
         int rows = this->mLDMatrices.find(key)->second.at(matIdx).rows();
         int cols = eqField.first.cols();

         // Create pointer to sparse operator
         const SparseMatrix * op = &this->mLDMatrices.find(key)->second.at(matIdx);

         // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
         eqField.first.block(eqStart, 0, rows, cols) += (*op)*linField.first.block(linStart, 0, rows, cols);
         eqField.second.block(eqStart, 0, rows, cols) += (*op)*linField.second.block(linStart, 0, rows, cols);
      }
   }

   void IPrognosticEquation::computeLinear(FieldComponents::Spectral::Id compId, DecoupledZMatrix& eqField, const int eqStart, SpectralFieldId fieldId, const MatrixZ& linField, const int linStart, const int matIdx) const
   {
      // Safety asserts
      assert(eqField.first.cols() == linField.cols());
      assert(eqField.second.cols() == linField.cols());

      // Create key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(compId, fieldId);

      if(this->mLZMatrices.count(key) > 0)
      {
         // Get number of rows and columns of block
         int rows = this->mLZMatrices.find(key)->second.at(matIdx).rows();
         int cols = eqField.first.cols();

         // Create pointer to sparse operator
         const SparseMatrixZ * op = &this->mLZMatrices.find(key)->second.at(matIdx);

         // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
         eqField.first.block(eqStart, 0, rows, cols) += op->real()*linField.block(linStart, 0, rows,cols).real();
         eqField.first.block(eqStart, 0, rows, cols) -= op->imag()*linField.block(linStart, 0, rows, cols).imag();
         eqField.second.block(eqStart, 0, rows, cols) += op->real()*linField.block(linStart, 0, rows, cols).imag();
         eqField.second.block(eqStart, 0, rows, cols) += op->imag()*linField.block(linStart, 0, rows, cols).real();
      }

      if(this->mLDMatrices.count(key) > 0)
      {
         // Get number of rows and columns of block
         int rows = this->mLDMatrices.find(key)->second.at(matIdx).rows();
         int cols = eqField.first.cols();

         // Create pointer to sparse operator
         const SparseMatrix * op = &this->mLDMatrices.find(key)->second.at(matIdx);

         // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
         eqField.first.block(eqStart, 0, rows, cols) += (*op)*linField.block(linStart, 0, rows, cols).real();
         eqField.second.block(eqStart, 0, rows, cols) += (*op)*linField.block(linStart, 0, rows, cols).imag();
      }
   }

   void IPrognosticEquation::computeLinear(FieldComponents::Spectral::Id compId, MatrixZ& eqField, const int eqStart, SpectralFieldId fieldId, const DecoupledZMatrix& linField, const int linStart, const int matIdx) const
   {
      // Safety asserts
      assert(eqField.cols() == linField.first.cols());
      assert(eqField.cols() == linField.second.cols());

      // Create key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(compId, fieldId);

      if(this->mLZMatrices.count(key) > 0)
      {
         // Get number of rows and columns of block
         int rows = this->mLZMatrices.find(key)->second.at(matIdx).rows();
         int cols = eqField.cols();

         // Create pointer to sparse operator
         const SparseMatrixZ * op = &this->mLZMatrices.find(key)->second.at(matIdx);

         // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
         eqField.block(eqStart, 0, rows, cols).real() += op->real()*linField.first.block(linStart, 0, rows, cols);
         eqField.block(eqStart, 0, rows, cols).real() -= op->imag()*linField.second.block(linStart, 0, rows, cols);
         eqField.block(eqStart, 0, rows, cols).imag() += op->imag()*linField.second.block(linStart, 0, rows, cols);
         eqField.block(eqStart, 0, rows, cols).imag() += op->imag()*linField.first.block(linStart, 0, rows, cols);
      }

      if(this->mLDMatrices.count(key) > 0)
      {
         // Get number of rows and columns of block
         int rows = this->mLDMatrices.find(key)->second.at(matIdx).rows();
         int cols = eqField.cols();

         // Create pointer to sparse operator
         const SparseMatrix * op = &this->mLDMatrices.find(key)->second.at(matIdx);

         // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
         eqField.block(eqStart, 0, rows, cols).real() += (*op)*linField.first.block(linStart, 0, rows, cols);
         eqField.block(eqStart, 0, rows, cols).imag() += (*op)*linField.second.block(linStart, 0, rows, cols);
      }
   }

   void IPrognosticEquation::computeLinear(FieldComponents::Spectral::Id compId, MatrixZ& eqField, const int eqStart, SpectralFieldId fieldId, const MatrixZ& linField, const int linStart, const int matIdx) const
   {
      // Safety asserts
      assert(eqField.cols() == linField.cols());

      // Create key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(compId, fieldId);

      if(this->mLZMatrices.count(key) > 0)
      {
         // Get number of rows and columns of block
         int rows = this->mLZMatrices.find(key)->second.at(matIdx).rows();
         int cols = eqField.cols();

         // Create pointer to sparse operator
         const SparseMatrixZ * op = &this->mLZMatrices.find(key)->second.at(matIdx);

         // Apply operator to field: eq += op*lin
         eqField.block(eqStart, 0, rows, cols) += (*op)*linField.block(linStart, 0, rows, cols);
      }

      if(this->mLDMatrices.count(key) > 0)
      {
         // Get number of rows and columns of block
         int rows = this->mLDMatrices.find(key)->second.at(matIdx).rows();
         int cols = eqField.cols();

         // Create pointer to sparse operator
         const SparseMatrix * op = &this->mLDMatrices.find(key)->second.at(matIdx);

         eqField.block(eqStart, 0, rows, cols).real() += (*op)*linField.block(linStart, 0, rows, cols).real();
         eqField.block(eqStart, 0, rows, cols).imag() += (*op)*linField.block(linStart, 0, rows, cols).imag();
      }
   }

   void IPrognosticEquation::applyNLQuasiInverse(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Get iterator to set of quasi-inverse matrices
      std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::const_iterator qIt = this->mNLMatrices.find(id);

      // Get number of rows
      int rows = this->couplingInfo(id).blockN(matIdx);

      // Safety asserts
      assert(qIt->second.at(matIdx).rows() == qIt->second.at(matIdx).cols());
      assert(qIt->second.at(matIdx).cols() == rows);

      // Multiply nonlinear term by quasi-inverse
      storage.first.block(start, 0, rows, storage.first.cols()) = qIt->second.at(matIdx)*storage.first.block(start, 0, rows, storage.first.cols());
      storage.second.block(start, 0, rows, storage.second.cols()) = qIt->second.at(matIdx)*storage.second.block(start, 0, rows, storage.second.cols());
   }

   void IPrognosticEquation::applyNLQuasiInverse(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start)
   {
      // Get iterator to set of quasi-inverse matrices
      std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::const_iterator qIt = this->mNLMatrices.find(id);

      // Get number of rows
      int rows = this->couplingInfo(id).blockN(matIdx);

      // Safety asserts
      assert(qIt->second.at(matIdx).rows() == qIt->second.at(matIdx).cols());
      assert(qIt->second.at(matIdx).cols() == rows);

      // Multiply nonlinear term by quasi-inverse
      storage.block(start, 0, rows, storage.cols()) = qIt->second.at(matIdx)*storage.block(start, 0, rows, storage.cols());
   }
}
}
