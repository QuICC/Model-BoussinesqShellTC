/** \file IEvolutionEquation.cpp
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
#include "Equations/IEvolutionEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IEvolutionEquation::IEvolutionEquation(SharedIEquationParameters spEqParams)
      : EquationData(spEqParams)
   {
   }

   IEvolutionEquation::~IEvolutionEquation()
   {
   }

   void IEvolutionEquation::init()
   {
      this->setCoupling();
   }

   void IEvolutionEquation::computeLinear(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const DecoupledZMatrix& linearField, const int matIdx, const int start) const
   {
//      // Get iterator to set of linear matrices
//      std::map<FieldComponents::Spectral::Id, std::vector<DecoupledZSparse> >::const_iterator lIt = this->mLMatrices.find(id);
//
//      // Get number of rows
//      int rows = this->couplingInfo(id).blockN(matIdx);
//
//      // Safety asserts
//      assert(lIt->second.at(matIdx).first.rows() == lIt->second.at(matIdx).first.cols());
//      assert(lIt->second.at(matIdx).second.rows() == lIt->second.at(matIdx).second.cols());
//      assert(lIt->second.at(matIdx).first.cols() == rows);
//      assert(lIt->second.at(matIdx).second.cols() == rows);
//
//      // Loop over all explicit linear terms
//      for(int k = 0; k < this->couplingInfo(FieldComponents::Spectral::SCALAR).nExplicit(); ++k)
//      {
//         storage.first.block(start, 0, rows, storage.first.cols()) += lIt->second.at(matIdx)*linearField.first.block(start, 0, rows, storage.first.cols());
//         storage.second.block(start, 0, rows, storage.second.cols()) += lIt->second.at(matIdx)*linearField.second.block(start, 0, rows, storage.second.cols());
//      }
   }

   void IEvolutionEquation::computeLinear(FieldComponents::Spectral::Id id, MatrixZ& storage, const MatrixZ& linearField, const int matIdx, const int start) const
   {
//      // Get iterator to set of linear matrices
//      std::map<FieldComponents::Spectral::Id, std::vector<DecoupledZSparse> >::const_iterator lIt = this->mLMatrices.find(id);
//
//      // Get number of rows
//      int rows = this->couplingInfo(id).blockN(matIdx);
//
//      // Safety asserts
//      assert(lIt->second.at(matIdx).first.rows() == lIt->second.at(matIdx).first.cols());
//      assert(lIt->second.at(matIdx).second.rows() == lIt->second.at(matIdx).second.cols());
//      assert(lIt->second.at(matIdx).first.cols() == rows);
//      assert(lIt->second.at(matIdx).second.cols() == rows);
//
//      // Loop over all explicit linear terms
//      CouplingInformation::field_iterator fIt;
//      CouplingInformation::field_iterator_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange();
//      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
//      {
//         storage.block(start, 0, rows, storage.cols()) += lIt->second.at(matIdx)*linearField.block(start, 0, rows, storage.cols());
//      }
   }

   void IEvolutionEquation::applyNLQuasiInverse(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start)
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

   void IEvolutionEquation::applyNLQuasiInverse(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start)
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
