/** \file IVectorEquation.cpp
 *  \brief Source of the base implementation of a vector equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IVectorEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IVectorEquation::IVectorEquation(SharedEquationParameters spEqParams)
      : IEquation(spEqParams)
   {
   }

   IVectorEquation::~IVectorEquation()
   {
   }

   void IVectorEquation::setUnknown(Datatypes::SharedVectorVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::VectorVariableType& IVectorEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::VectorVariableType& IVectorEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   void IVectorEquation::updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId)
   {
      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().cols() < rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().rComp(compId).setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().rows()));
   }

   void IVectorEquation::storeSolution(FieldComponents::Spectral::Id compId, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      if(this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

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
               this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(MHDComplex(storage.first(k),storage.second(k)),i,j,matIdx);

               // increase storage counter
               k++;
            }
         }
      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         // Get mode indexes
         ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = this->unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(MHDComplex(storage.first(k),storage.second(k)),i,mode(1),mode(0));

            // increase storage counter
            k++;
         }
      }
   }

   void IVectorEquation::storeSolution(FieldComponents::Spectral::Id compId, const MatrixZ& storage, const int matIdx, const int start)
   {
      if(this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

         //Safety assertion
         assert(rows*cols+start <= storage.rows());

         // Copy data
         int k = start;
         for(int j = 0; j < cols; j++)
         {
            for(int i = 0; i < rows; i++)
            {
               // Copy timestep output into field
               this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(storage(k),i,j,matIdx);

               // increase storage counter
               k++;
            }
         }
      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         // Get mode indexes
         ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = this->unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(storage(k),i,mode(1),mode(0));

            // increase storage counter
            k++;
         }
      }
   }

   void copyUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // matIdx is the index of the slowest varying direction
      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

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
               storage.first(k) = eq.unknown().dom(0).perturbation().comp(compId).point(i,j,matIdx).real();

               // Copy field imaginary value into storage
               storage.second(k) = eq.unknown().dom(0).perturbation().comp(compId).point(i,j,matIdx).imag();

               // increase storage counter
               k++;
            }
         }

      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         //Safety assertion
         assert(start >= 0);
         assert(start < storage.first.size());
         assert(start < storage.second.size());

         // Get mode indexes
         ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy field real value into storage
            storage.first(k) = eq.unknown().dom(0).perturbation().comp(compId).point(i,mode(1),mode(0)).real();

            // Copy field imaginary value into storage
            storage.second(k) = eq.unknown().dom(0).perturbation().comp(compId).point(i,mode(1),mode(0)).imag();

            // increase storage counter
            k++;
         }
      }
   }

   void copyUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start)
   {
      // matIdx is the index of the slowest varying direction
      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

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
               storage(k) = eq.unknown().dom(0).perturbation().comp(compId).point(i,j,matIdx);

               // increase storage counter
               k++;
            }
         }

      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         //Safety assertion
         assert(start >= 0);
         assert(start < storage.size());

         // Get mode indexes
         ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy field value into storage
            storage(k) = eq.unknown().dom(0).perturbation().comp(compId).point(i,mode(1),mode(0));

            // increase storage counter
            k++;
         }
      }
   }

   void copyNonlinear(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Check if a nonlinear computation took place
      if(eq.couplingInfo(compId).hasNonlinear())
      {
         // simply copy values from unknown
         copyUnknown(eq, compId, storage, matIdx, start);

      // Without nonlinear computation the values have to be initialised to zero   
      } else
      {
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

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
                  storage.first(k) = 0.0;

                  // Copy field imaginary value into storage
                  storage.second(k) = 0.0;

                  // increase storage counter
                  k++;
               }
            }
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            //Safety assertion
            assert(start >= 0);
            assert(start < storage.first.size());
            assert(start < storage.second.size());

            // Get mode indexes
            ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

            // Copy data
            int k = start;
            for(int i = 0; i < rows; i++)
            {
               // Copy field real value into storage
               storage.first(k) = 0.0;

               // Copy field imaginary value into storage
               storage.second(k) = 0.0;

               // increase storage counter
               k++;
            }
         }
      }
   }

   void copyNonlinear(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start)
   {
      // Check if a nonlinear computation took place
      if(eq.couplingInfo(compId).hasNonlinear())
      {
         // simply copy values from unknown
         copyUnknown(eq, compId, storage, matIdx, start);

      // Without nonlinear computation the values have to be initialised to zero   
      } else
      {
         // matIdx is the index of the slowest varying direction
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

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
                  storage(k) = 0.0;

                  // increase storage counter
                  k++;
               }
            }

         // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            //Safety assertion
            assert(start >= 0);
            assert(start < storage.size());

            // Get mode indexes
            ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

            // Copy data
            int k = start;
            for(int i = 0; i < rows; i++)
            {
               // Copy field value into storage
               storage(k) = 0.0;

               // increase storage counter
               k++;
            }
         }
      }
   }

   void addSource(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Add source term if required
      if(eq.couplingInfo(compId).hasSource())
      {
         // matIdx is the index of the slowest varying direction
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

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
                  // Get source value
                  MHDComplex source = eq.sourceTerm(compId, i, j, matIdx);

                  // Add real part of source term
                  storage.first(k) += source.real();

                  // Add imaginary part of source term
                  storage.second(k) += source.imag();

                  // increase storage counter
                  k++;
               }
            }

         // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            //Safety assertion
            assert(start >= 0);
            assert(start < storage.first.size());
            assert(start < storage.second.size());

            // Get mode indexes
            ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

            // Copy data
            int k = start;
            for(int i = 0; i < rows; i++)
            {
               // Get source value
               MHDComplex source = eq.sourceTerm(compId, i, mode(1), mode(0));

               // Add real part of source term
               storage.first(k) += source.real();

               // Add imaginary part of source term
               storage.second(k) += source.imag();

               // increase storage counter
               k++;
            }
         }
      }
   }

   void addSource(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start)
   {
      // Add source term if required
      if(eq.couplingInfo(compId).hasSource())
      {
         // matIdx is the index of the slowest varying direction
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

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
                  // Add source term
                  storage(k) = eq.sourceTerm(compId, i, j, matIdx);

                  // increase storage counter
                  k++;
               }
            }

         // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            //Safety assertion
            assert(start >= 0);
            assert(start < storage.size());

            // Get mode indexes
            ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

            // Copy data
            int k = start;
            for(int i = 0; i < rows; i++)
            {
               // Add source term
               storage(k) = eq.sourceTerm(compId, i, mode(1), mode(0));

               // increase storage counter
               k++;
            }
         }
      }
   }
}
}
