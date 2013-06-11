/** \file IScalarPEquation.cpp
 *  \brief Source of the base implementation of a scalar prognostic equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IScalarPEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IScalarPEquation::IScalarPEquation(SharedEquationParameters spEqParams)
      : IPrognosticEquation(spEqParams)
   {
   }

   IScalarPEquation::~IScalarPEquation()
   {
   }

   void IScalarPEquation::setUnknown(Datatypes::SharedScalarVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::ScalarVariableType& IScalarPEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::ScalarVariableType& IScalarPEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   void IScalarPEquation::prepareTimestep(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id id)
   {
      // Assert dealiasing has taken place!
      assert(id == FieldComponents::Spectral::SCALAR);
      assert(this->rUnknown().rDom(0).rPerturbation().data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().data().cols() == rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));
   }

   void IScalarPEquation::timestepInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Copy field values into timestep input
      this->copyTInput(id, storage, matIdx, start);

      // Apply quasi-inverse to nonlinear terms
      this->applyNLQuasiInverse(id, storage, matIdx, start);
   }

   void IScalarPEquation::timestepInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start)
   {
      // Copy field values into timestep input
      this->copyTInput(id, storage, matIdx, start);

      // Apply quasi-inverse to nonlinear terms
      this->applyNLQuasiInverse(id, storage, matIdx, start);
   }

   void IScalarPEquation::timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Copy timestep output into field value
      this->copyTOutput(id, storage, matIdx, start);
   }

   void IScalarPEquation::timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      // Copy timestep output into field value
      this->copyTOutput(id, storage, matIdx, start);
   }

   void IScalarPEquation::copyTInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      if(this->couplingInfo(id).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);
         assert(start < storage.first.size());
         assert(start < storage.second.size());
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
      } else if(this->couplingInfo(id).indexType() == CouplingInformation::MODE)
      {
         //Safety assertion
         assert(start >= 0);
         assert(start < storage.first.size());
         assert(start < storage.second.size());

         // Get mode indexes
         ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = this->unknown().dom(0).perturbation().slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy field real value into storage
            storage.first(k) = this->unknown().dom(0).perturbation().point(i,mode(1),mode(0)).real();

            // Copy field imaginary value into storage
            storage.second(k) = this->unknown().dom(0).perturbation().point(i,mode(1),mode(0)).imag();

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarPEquation::copyTInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start)
   {
      if(this->couplingInfo(id).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);
         assert(start < storage.size());
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
      } else if(this->couplingInfo(id).indexType() == CouplingInformation::MODE)
      {
         //Safety assertion
         assert(start >= 0);
         assert(start < storage.size());

         // Get mode indexes
         ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

         int rows = this->unknown().dom(0).perturbation().slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy field value into storage
            storage(k) = this->unknown().dom(0).perturbation().point(i,mode(1),mode(0));

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarPEquation::copyTOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      if(this->couplingInfo(id).indexType() == CouplingInformation::SLOWEST)
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
      } else if(this->couplingInfo(id).indexType() == CouplingInformation::MODE)
      {
         // Get mode indexes
         ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = this->unknown().dom(0).perturbation().slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            this->rUnknown().rDom(0).rPerturbation().setPoint(MHDComplex(storage.first(k),storage.second(k)),i,mode(1),mode(0));

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarPEquation::copyTOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      if(this->couplingInfo(id).indexType() == CouplingInformation::SLOWEST)
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
      } else if(this->couplingInfo(id).indexType() == CouplingInformation::MODE)
      {
         // Get mode indexes
         ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = this->unknown().dom(0).perturbation().slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            this->rUnknown().rDom(0).rPerturbation().setPoint(storage(k),i,mode(1),mode(0));

            // increase storage counter
            k++;
         }
      }
   }
}
}
