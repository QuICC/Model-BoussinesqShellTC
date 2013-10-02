/** 
 * @file IScalarEquation.cpp
 * @brief Source of the base implementation of a scalar equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Base/MathConstants.hpp"
#include "TypeSelectors/EquationToolsSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IScalarEquation::IScalarEquation(SharedEquationParameters spEqParams)
      : IEquation(spEqParams)
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

   SharedResolution IScalarEquation::spRes() const
   {
      return this->unknown().dom(0).spRes();
   }

   void IScalarEquation::updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);
      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().data().cols() == rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));
   }

   void IScalarEquation::storeSolution(FieldComponents::Spectral::Id compId, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      if(this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

         //Safety assertion
         assert(rows*cols+start <= storage.real().rows());
         assert(rows*cols+start <= storage.imag().rows());

         // Copy data
         int k = start;
         for(int j = 0; j < cols; j++)
         {
            for(int i = 0; i < rows; i++)
            {
               // Copy timestep output into field
               this->rUnknown().rDom(0).rPerturbation().setPoint(MHDComplex(storage.real()(k),storage.imag()(k)),i,j,matIdx);

               // increase storage counter
               k++;
            }
         }
      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         // Get mode indexes
         ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = this->unknown().dom(0).perturbation().slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            this->rUnknown().rDom(0).rPerturbation().setPoint(MHDComplex(storage.real()(k),storage.imag()(k)),i,mode(1),mode(0));

            // increase storage counter
            k++;
         }
      }
   }

   void IScalarEquation::storeSolution(FieldComponents::Spectral::Id compId, const MatrixZ& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      if(this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
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
      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::MODE)
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

   void IScalarEquation::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Make sure it is safe to do nothing
      bool needInit = this->couplingInfo(FieldComponents::Spectral::SCALAR).hasQuasiInverse();

      CouplingInformation::FieldId_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange();
      needInit = needInit && (fRange.first == fRange.second);

      // Initialise spectral matrices
      if(needInit)
      {
         if(EquationToolsType::EIGEN_DIMS == 0)
         {
            this->initSpectralMatricesNoEigen(spBcIds, FieldComponents::Spectral::SCALAR);
         } else if(EquationToolsType::EIGEN_DIMS == 1)
         {
            this->initSpectralMatrices1DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
         } else if(EquationToolsType::EIGEN_DIMS == 2)
         {
            this->initSpectralMatrices2DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
         } else if(EquationToolsType::EIGEN_DIMS == 3)
         {
            this->initSpectralMatrices3DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
         } else
         {
            throw Exception("initSpectralMatrices: unhandled number of eigen dimensions!");
         }
      }
   }

   void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      // matIdx is the index of the slowest varying direction
      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);
         assert(start < storage.real().size());
         assert(start < storage.imag().size());
         assert(rows*cols+start <= storage.real().rows());
         assert(rows*cols+start <= storage.imag().rows());

         // Copy data
         int k = start;
         for(int j = 0; j < cols; j++)
         {
            for(int i = 0; i < rows; i++)
            {
               // Copy field real value into storage
               storage.real()(k) = eq.unknown().dom(0).perturbation().point(i,j,matIdx).real();

               // Copy field imaginary value into storage
               storage.imag()(k) = eq.unknown().dom(0).perturbation().point(i,j,matIdx).imag();

               // increase storage counter
               k++;
            }
         }

      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         //Safety assertion
         assert(start >= 0);
         assert(start < storage.real().size());
         assert(start < storage.imag().size());

         // Get mode indexes
         ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = eq.unknown().dom(0).perturbation().slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy field real value into storage
            storage.real()(k) = eq.unknown().dom(0).perturbation().point(i,mode(1),mode(0)).real();

            // Copy field imaginary value into storage
            storage.imag()(k) = eq.unknown().dom(0).perturbation().point(i,mode(1),mode(0)).imag();

            // increase storage counter
            k++;
         }
      }
   }

   void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      // matIdx is the index of the slowest varying direction
      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();

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
               storage(k) = eq.unknown().dom(0).perturbation().point(i,j,matIdx);

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

         int rows = eq.unknown().dom(0).perturbation().slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = 0; i < rows; i++)
         {
            // Copy field value into storage
            storage(k) = eq.unknown().dom(0).perturbation().point(i,mode(1),mode(0));

            // increase storage counter
            k++;
         }
      }
   }

   void copyNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

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
            int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();

            //Safety assertion
            assert(start >= 0);
            assert(start < storage.real().size());
            assert(start < storage.imag().size());
            assert(rows*cols+start <= storage.real().rows());
            assert(rows*cols+start <= storage.imag().rows());

            // Set data to zero
            int k = start;
            for(int j = 0; j < cols; j++)
            {
               for(int i = 0; i < rows; i++)
               {
                  // Copy field real value into storage
                  storage.real()(k) = 0.0;

                  // Copy field imaginary value into storage
                  storage.imag()(k) = 0.0;

                  // increase storage counter
                  k++;
               }
            }

            // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            //Safety assertion
            assert(start >= 0);
            assert(start < storage.real().size());
            assert(start < storage.imag().size());

            // Get mode indexes
            ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
            int rows = eq.unknown().dom(0).perturbation().slice(mode(0)).rows();

            // Set data to zero
            int k = start;
            for(int i = 0; i < rows; i++)
            {
               // Copy field real value into storage
               storage.real()(k) = 0.0;

               // Copy field imaginary value into storage
               storage.imag()(k) = 0.0;

               // increase storage counter
               k++;
            }
         }
      }
   }

   void copyNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

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
            int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();

            //Safety assertion
            assert(start >= 0);
            assert(start < storage.size());
            assert(rows*cols+start <= storage.rows());

            // Set data to zero
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

            int rows = eq.unknown().dom(0).perturbation().slice(mode(0)).rows();

            // Set data to zero
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

   void addSource(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      // Add source term if required
      if(eq.couplingInfo(compId).hasSource())
      {
         // matIdx is the index of the slowest varying direction
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();

            //Safety assertion
            assert(start >= 0);
            assert(start < storage.real().size());
            assert(start < storage.imag().size());
            assert(rows*cols+start <= storage.real().rows());
            assert(rows*cols+start <= storage.imag().rows());

            // Copy data
            int k = start;
            for(int j = 0; j < cols; j++)
            {
               for(int i = 0; i < rows; i++)
               {
                  // Get source value
                  MHDComplex source = eq.sourceTerm(compId, i, j, matIdx);

                  // Add real part of source term
                  storage.real()(k) += source.real();

                  // Add imaginary part of source term
                  storage.imag()(k) += source.imag();

                  // increase storage counter
                  k++;
               }
            }

         // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            //Safety assertion
            assert(start >= 0);
            assert(start < storage.real().size());
            assert(start < storage.imag().size());

            // Get mode indexes
            ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
            int rows = eq.unknown().dom(0).perturbation().slice(mode(0)).rows();

            // Copy data
            int k = start;
            for(int i = 0; i < rows; i++)
            {
               // Get source value
               MHDComplex source = eq.sourceTerm(compId, i, mode(1), mode(0));

               // Add real part of source term
               storage.real()(k) += source.real();

               // Add imaginary part of source term
               storage.imag()(k) += source.imag();

               // increase storage counter
               k++;
            }
         }
      }
   }

   void addSource(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      // Add source term if required
      if(eq.couplingInfo(compId).hasSource())
      {
         // matIdx is the index of the slowest varying direction
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();

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

            int rows = eq.unknown().dom(0).perturbation().slice(mode(0)).rows();

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
