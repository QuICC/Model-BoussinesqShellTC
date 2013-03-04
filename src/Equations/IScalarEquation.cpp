/** \file IScalarEquation.cpp
 *  \brief Source of the base implementation of a scalar equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IScalarEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IScalarEquation::IScalarEquation(SharedIEquationParameters spEqParams)
      : IEvolutionEquation(spEqParams)
   {
   }

   IScalarEquation::~IScalarEquation()
   {
   }

   void IScalarEquation::setUnknown(Datatypes::SharedScalarVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::ScalarVariableType& IScalarEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::ScalarVariableType& IScalarEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   void IScalarEquation::computeLinear(Datatypes::SpectralScalarType& rRHS) const
   {
      // Empty default implementation
   }

   void IScalarEquation::prepareTimestep(const Datatypes::SpectralScalarType& rhs)
   {
      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().data().cols() < rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));
   }

   void IScalarEquation::timestepInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Copy field values into timestep input
      this->copyTInput(id,storage,matIdx,start);

      // Apply quasi-inverse
      this->applyQuasiInverse(id,storage,matIdx,start);
   }

   void IScalarEquation::timestepInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start)
   {
      // Copy field values into timestep input
      this->copyTInput(id,storage,matIdx,start);

      // Apply quasi-inverse
      this->applyQuasiInverse(id,storage,matIdx,start);
   }

   void IScalarEquation::timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Copy timestep output into field value
      this->copyTOutput(id, storage, matIdx, start);
   }

   void IScalarEquation::timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      // Copy timestep output into field value
      this->copyTOutput(id, storage, matIdx, start);
   }

   void IScalarEquation::applyQuasiInverse(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
      int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();
      int inCols = storage.first.cols();

      storage.first.block(start,0,rows*cols,inCols) = this->mNLMatrices.find(id)->second.at(matIdx)*storage.first.block(start,0,rows*cols,inCols);
      storage.second.block(start,0,rows*cols,inCols) = this->mNLMatrices.find(id)->second.at(matIdx)*storage.second.block(start,0,rows*cols,inCols);
   }

   void IScalarEquation::applyQuasiInverse(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start)
   {
      int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
      int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();
      int inCols = storage.cols();

      storage.block(start,0,rows*cols,inCols) = this->mNLMatrices.find(id)->second.at(matIdx)*storage.block(start,0,rows*cols,inCols);
   }

   void IScalarEquation::copyTInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
      int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

      //Safety assertion
      assert(rows*cols+start <= storage.first.rows());
      assert(rows*cols+start <= storage.second.rows());

      // Copy data
      int k = start;
      for(int j = 0; j < cols; j++)
      {
         for(int i = 0; i < rows; i++)
         {
            // Copy field real value into storage
            storage.first(k) = this->unknown().dom(0).perturbation().point(i,j,matIdx).real();

            // Copy field imaginary value into storage
            storage.second(k) = this->unknown().dom(0).perturbation().point(i,j,matIdx).imag();

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarEquation::copyTInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start)
   {
      int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
      int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

      //Safety assertion
      assert(rows*cols+start <= storage.rows());

      // Copy data
      int k = start;
      for(int j = 0; j < cols; j++)
      {
         for(int i = 0; i < rows; i++)
         {
            // Copy field into storage
            storage(k) = this->unknown().dom(0).perturbation().point(i,j,matIdx);

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarEquation::copyTOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
      int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

      //Safety assertion
      assert(rows*cols+start <= storage.first.rows());
      assert(rows*cols+start <= storage.second.rows());

      // Copy data
      int k = start;
      for(int j = 0; j < cols; j++)
      {
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            this->rUnknown().rDom(0).rPerturbation().setPoint(MHDComplex(storage.first(k),storage.second(k)),i,j,matIdx);

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarEquation::copyTOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
      int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

      //Safety assertion
      assert(rows*cols+start <= storage.rows());

      // Copy data
      int k = start;
      for(int j = 0; j < cols; j++)
      {
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            this->rUnknown().rDom(0).rPerturbation().setPoint(storage(k),i,j,matIdx);

            // increase storage counter
            k++;
         }
      }
   }
}
}
