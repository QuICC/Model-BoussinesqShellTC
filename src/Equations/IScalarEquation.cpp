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

   IScalarEquation::IScalarEquation(SharedEquationParameters spEqParams)
      : IEvolutionEquation(spEqParams)
   {
   }

   void IScalarEquation::applyQuasiInverse(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      int rows = 2*this->unknown().dom(0).perturbation().slice(matIdx).rows()/3;
      int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();
      int inCols = storage.first.cols();

      storage.first.block(start,0,rows*cols,inCols) = this->mNLMatrices.find(id)->second.at(matIdx)*storage.first.block(start,0,rows*cols,inCols);
      storage.second.block(start,0,rows*cols,inCols) = this->mNLMatrices.find(id)->second.at(matIdx)*storage.second.block(start,0,rows*cols,inCols);
   }

   void IScalarEquation::applyQuasiInverse(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start)
   {
      int rows = 2*this->unknown().dom(0).perturbation().slice(matIdx).rows()/3;
      int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();
      int inCols = storage.cols();

      storage.block(start,0,rows*cols,inCols) = this->mNLMatrices.find(id)->second.at(matIdx)*storage.block(start,0,rows*cols,inCols);
   }

   void IScalarEquation::copyTInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      int rows = 2*this->unknown().dom(0).perturbation().slice(matIdx).rows()/3;
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
            storage.first(k) = this->unknown().dom(0).perturbation().slice(matIdx)(i,j).real();

            // Copy field imaginary value into storage
            storage.second(k) = this->unknown().dom(0).perturbation().slice(matIdx)(i,j).imag();

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarEquation::copyTInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start)
   {
      int rows = 2*this->unknown().dom(0).perturbation().slice(matIdx).rows()/3;
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
            storage(k) = this->unknown().dom(0).perturbation().slice(matIdx)(i,j);

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarEquation::copyTOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      int rows = 2*this->unknown().dom(0).perturbation().slice(matIdx).rows()/3;
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
            // Copy timestep output into field real component
            this->rUnknown().rDom(0).rPerturbation().rSlice(matIdx)(i,j).real() = storage.first(k);

            // Copy timestep output into field imaginary component
            this->rUnknown().rDom(0).rPerturbation().rSlice(matIdx)(i,j).imag() = storage.second(k);

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarEquation::copyTOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      int rows = 2*this->unknown().dom(0).perturbation().slice(matIdx).rows()/3;
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
            this->rUnknown().rDom(0).rPerturbation().rSlice(matIdx)(i,j) = storage(k);

            // increase storage counter
            k++;
         }
      }
//      std::cerr << rows << " x " << cols << std::endl;
//      std::cerr << " --------------------------------------------------------------------- " << std::endl;
//      std::cerr << storage.transpose() << std::endl;
//      std::cerr << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << std::endl;
//      std::cerr << this->rUnknown().rDom(0).rPerturbation().rSlice(matIdx) << std::endl;
//      std::cerr << " --------------------------------------------------------------------- " << std::endl;
   }
}
