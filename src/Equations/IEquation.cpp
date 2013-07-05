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
#include "Base/MathConstants.hpp"

#include <iostream>
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

   void IEquation::initSpectralMatrices1DEigen(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Get the number of systems
      int nSystems = this->couplingInfo(compId).nSystems();

      // Boxscale
      MHDFloat boxScale = this->spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);

      //
      // Initialise the quasi-inverse operators for the nonlinear terms (if required)
      //
      if(this->couplingInfo(compId).hasQuasiInverse())
      {
         this->mNLMatrices.insert(std::make_pair(compId, std::vector<SparseMatrix>()));
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::iterator qIt = this->mNLMatrices.find(compId);
         qIt->second.reserve(nSystems);
         for(int i = 0; i < nSystems; ++i)
         {
            qIt->second.push_back(SparseMatrix());

            this->setQuasiInverse(compId, qIt->second.back());
         }
      }

      //
      // Initialise the explicit linear operators
      //
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = this->couplingInfo(compId).explicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         std::vector<DecoupledZSparse> tmpMat;
         tmpMat.reserve(nSystems);

         bool isComplex = false;

         // Create matrices
         for(int i = 0; i < nSystems; ++i)
         {
            MHDFloat k_ = boxScale*static_cast<MHDFloat>(this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i));

            // Get linear block
            tmpMat.push_back(DecoupledZSparse());
            this->setExplicitLinearBlock(compId, tmpMat.at(i), *fIt, k_);

            // Explicit operator requires an additional minus sign
            tmpMat.at(i).first = -tmpMat.at(i).first;
            tmpMat.at(i).second = -tmpMat.at(i).second;

            isComplex = isComplex || (tmpMat.at(i).second.nonZeros() > 0);
         }

         // Create key
         std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(compId, *fIt);

         // Select real or complex operator
         if(isComplex)
         {
            this->mLZMatrices.insert(std::make_pair(key, std::vector<SparseMatrixZ>()));
            this->mLZMatrices.find(key)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               SparseMatrixZ tmp = tmpMat.at(i).first.cast<MHDComplex>() + MathConstants::cI*tmpMat.at(i).second;
               this->mLZMatrices.find(key)->second.push_back(tmp);
            }
         } else
         {
            this->mLDMatrices.insert(std::make_pair(key, std::vector<SparseMatrix>()));
            this->mLDMatrices.find(key)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               this->mLDMatrices.find(key)->second.push_back(SparseMatrix());

               this->mLDMatrices.find(key)->second.back().swap(tmpMat.at(i).first);
            }
         }
      }
   }

   void IEquation::initSpectralMatrices2DEigen(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Get the number of systems
      int nSystems = this->couplingInfo(compId).nSystems();

      // Boxscale
      MHDFloat box2DScale = this->spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);
      MHDFloat box3DScale = this->spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D);

      //
      // Initialise the quasi-inverse operators for the nonlinear terms (if required)
      //
      if(this->couplingInfo(compId).hasQuasiInverse())
      {
         this->mNLMatrices.insert(std::make_pair(compId, std::vector<SparseMatrix>()));
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::iterator qIt = this->mNLMatrices.find(compId);
         qIt->second.reserve(nSystems);
         for(int i = 0; i < nSystems; ++i)
         {
            qIt->second.push_back(SparseMatrix());

            this->setQuasiInverse(compId, qIt->second.back());
         }
      }

      //
      // Initialise the explicit linear operators
      //
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = this->couplingInfo(compId).explicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         std::vector<DecoupledZSparse> tmpMat;
         tmpMat.reserve(nSystems);

         bool isComplex = false;

         // Create matrices
         for(int i = 0; i < nSystems; ++i)
         {
            // Get mode indexes
            ArrayI mode = this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(i);

            // Get 2D wave number rescaled to box size
            MHDFloat k2D_ = box2DScale*static_cast<MHDFloat>(mode(0));
            // Get 3D wave number rescaled to box size
            MHDFloat k3D_ = box3DScale*static_cast<MHDFloat>(mode(1));

            // Get linear block
            tmpMat.push_back(DecoupledZSparse());
            this->setExplicitLinearBlock(compId, tmpMat.at(i), *fIt, k2D_, k3D_);

            // Explicit operator requires an additional minus sign
            tmpMat.at(i).first = -tmpMat.at(i).first;
            tmpMat.at(i).second = -tmpMat.at(i).second;

            isComplex = isComplex || (tmpMat.at(i).second.nonZeros() > 0);
         }

         // Create key
         std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(compId, *fIt);

         // Select real or complex operator
         if(isComplex)
         {
            this->mLZMatrices.insert(std::make_pair(key, std::vector<SparseMatrixZ>()));
            this->mLZMatrices.find(key)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               SparseMatrixZ tmp = tmpMat.at(i).first.cast<MHDComplex>() + MathConstants::cI*tmpMat.at(i).second;
               this->mLZMatrices.find(key)->second.push_back(tmp);
            }
         } else
         {
            this->mLDMatrices.insert(std::make_pair(key, std::vector<SparseMatrix>()));
            this->mLDMatrices.find(key)->second.back().reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               this->mLDMatrices.find(key)->second.push_back(SparseMatrix());

               this->mLDMatrices.find(key)->second.back().swap(tmpMat.at(i).first);
            }
         }
      }
   }

   void IEquation::setQuasiInverse(FieldComponents::Spectral::Id comp, SparseMatrix &mat) const
   {
      throw Exception("setQuasiInverse: dummy implementation was called!");
   }

   void IEquation::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      throw Exception("setExplicitLinearBlock: dummy implementation was called!");
   }

   void IEquation::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k2D, const MHDFloat k3D) const
   {
      throw Exception("setExplicitLinearBlock: dummy implementation was called!");
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

   void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx)
   {
      // Compute with complex linear operator
      if(eq.hasExplicitZLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.explicitZLinear(compId, fieldId, matIdx);

         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            // Get number of rows and columns of block
            int rows = op->rows()*op->cols();
            int cols = eqField.first.cols();

            //Safety assertion
            assert(rows+eqStart <= eqField.first.rows());
            assert(rows+eqStart <= eqField.second.rows());

            /// \mhdBug very bad and slow implementation!
            ArrayZ   tmp(rows);
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

            // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
            eqField.first.block(eqStart, 0, rows, cols) += op->real()*tmp.real();
            eqField.first.block(eqStart, 0, rows, cols) -= op->imag()*tmp.imag();
            eqField.second.block(eqStart, 0, rows, cols) += op->real()*tmp.imag();
            eqField.second.block(eqStart, 0, rows, cols) += op->imag()*tmp.real();

         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            // Get mode indexes
            ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

            // Assert correct sizes
            assert(op->rows() == explicitField.slice(mode(0)).rows());
            assert(eqField.first.cols() == 1);
            assert(eqField.second.cols() == 1);

            // Get number of rows and columns of block
            int rows = op->rows();
            int cols = eqField.first.cols();

            // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
            eqField.first.block(eqStart, 0, rows, cols) += op->real()*explicitField.slice(mode(0)).col(mode(1)).real();
            eqField.first.block(eqStart, 0, rows, cols) -= op->imag()*explicitField.slice(mode(0)).col(mode(1)).imag();
            eqField.second.block(eqStart, 0, rows, cols) += op->real()*explicitField.slice(mode(0)).col(mode(1)).imag();
            eqField.second.block(eqStart, 0, rows, cols) += op->imag()*explicitField.slice(mode(0)).col(mode(1)).real();
         }
      }

      // Compute with real linear operator
      if(eq.hasExplicitDLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.explicitDLinear(compId, fieldId, matIdx);

         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            // Get number of rows and columns of block
            int rows = explicitField.slice(matIdx).size();
            int cols = eqField.first.cols();

            //Safety assertion
            assert(rows+eqStart <= eqField.first.rows());
            assert(rows+eqStart <= eqField.second.rows());

            /// \mhdBug very bad and slow implementation!
            ArrayZ   tmp(rows);
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

            // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
            eqField.first.block(eqStart, 0, rows, cols) += (*op)*tmp.real();
            eqField.second.block(eqStart, 0, rows, cols) += (*op)*tmp.imag();

         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            // Get mode indexes
            ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

            // Assert correct sizes
            assert(op->rows() == explicitField.slice(mode(0)).rows());
            assert(eqField.first.cols() == 1);
            assert(eqField.second.cols() == 1);

            // Get number of rows and columns of block
            int rows = op->rows();
            int cols = eqField.first.cols();

            // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
            eqField.first.block(eqStart, 0, rows, cols) += (*op)*explicitField.slice(mode(0)).col(mode(1)).real();
            eqField.second.block(eqStart, 0, rows, cols) += (*op)*explicitField.slice(mode(0)).col(mode(1)).imag();
         }
      }
   }

   void addExplicitLinear(const IEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& eqField, const int eqStart, SpectralFieldId fieldId, const Datatypes::SpectralScalarType& explicitField, const int matIdx)
   {
      // Compute with complex linear operator
      if(eq.hasExplicitZLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.explicitZLinear(compId, fieldId, matIdx);

         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            // Get number of rows and columns of block
            int rows = op->rows()*op->cols();
            int cols = eqField.cols();

            //Safety assertion
            assert(rows+eqStart <= eqField.rows());

            /// \mhdBug very bad and slow implementation!
            ArrayZ   tmp(rows);
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

            // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
            eqField.block(eqStart, 0, rows, cols) += (*op)*tmp;

         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            // Get mode indexes
            ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

            // Assert correct sizes
            assert(op->rows() == explicitField.slice(mode(0)).rows());
            assert(eqField.cols() == 1);

            // Get number of rows and columns of block
            int rows = op->rows();
            int cols = eqField.cols();

            // Apply operator to field: Re(eq) += Re(op)*Re(lin) - Im(op)*Im(lin), Im(eq) += Re(op)*Im(lin) + Im(op)*Re(lin)
            eqField.block(eqStart, 0, rows, cols) += (*op)*explicitField.slice(mode(0)).col(mode(1));
         }
      }

      // Compute with real linear operator
      if(eq.hasExplicitDLinear(compId, fieldId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.explicitDLinear(compId, fieldId, matIdx);

         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            // Get number of rows and columns of block
            int rows = explicitField.slice(matIdx).size();
            int cols = eqField.cols();

            //Safety assertion
            assert(rows+eqStart <= eqField.rows());

            /// \mhdBug very bad and slow implementation!
            ArrayZ   tmp(rows);
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

            // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
            eqField.block(eqStart, 0, rows, cols).real() += (*op)*tmp.real();
            eqField.block(eqStart, 0, rows, cols).imag() += (*op)*tmp.imag();

         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            // Get mode indexes
            ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

            // Assert correct sizes
            assert(op->rows() == explicitField.slice(mode(0)).rows());
            assert(eqField.cols() == 1);

            // Get number of rows and columns of block
            int rows = op->rows();
            int cols = eqField.cols();

            // Apply operator to field: Re(eq) += op*Re(lin), Im(eq) += op*Im(lin)
            eqField.block(eqStart, 0, rows, cols).real() += (*op)*explicitField.slice(mode(0)).col(mode(1)).real();
            eqField.block(eqStart, 0, rows, cols).imag() += (*op)*explicitField.slice(mode(0)).col(mode(1)).imag();
         }
      }
   }

   void quasiInverseBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      throw Exception("quasiInverseBlock: dummy implementation called!");
   }

   void timeBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat k)
   {
      throw Exception("timeBlock: dummy implementation called!");
   }

   void linearBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      throw Exception("linearBlock: dummy implementation called!");
   }

   void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      throw Exception("boundaryBlock: dummy implementation called!");
   }
}
}
