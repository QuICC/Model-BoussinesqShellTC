/** \file IBoussinesqPerBetaCylGScalarEquation.cpp
 *  \brief Source of the implementation of the streamfunction equation in the 3DQG beta model with periodic radius
 */

// Configuration includes
//

// System includes
//

// External includes
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/IBoussinesqPerBetaCylGScalarEquation.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGSystem.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IBoussinesqPerBetaCylGScalarEquation::IBoussinesqPerBetaCylGScalarEquation(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   IBoussinesqPerBetaCylGScalarEquation::~IBoussinesqPerBetaCylGScalarEquation()
   {
   }

   void IBoussinesqPerBetaCylGScalarEquation::setCoupling()
   {
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Get Y dimension
      int nY = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();

      int modes = 0;
      for(int i = 0; i < nY; i++)
      {
         modes += this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      // Set coupling information
      this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);
      BoussinesqPerBetaCylGSystem::setCouplingInfo(this->mCouplingInfos.find(FieldComponents::Spectral::SCALAR)->second, eqId, nZ, modes);
   }

   DecoupledZSparse IBoussinesqPerBetaCylGScalarEquation::linearRow(FieldComponents::Spectral::Id comp, const int matIdx) const
   {
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Get mode indexes
      ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get radial wave number rescaled to box size
      MHDFloat kX_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1));
      // Get azimuthal wave number rescaled to box size
      MHDFloat kY_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0));

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = this->eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat Gamma = this->eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = this->eqParams().nd(NonDimensional::CHI);

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)), SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)),SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)));
      SparseMatrix  tmp(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx));

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::field_iterator fIt;
      CouplingInformation::field_iterator_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix   blockMatrix(this->couplingInfo(comp).nBlocks(),this->couplingInfo(comp).nBlocks());
         blockMatrix.insert(this->couplingInfo(comp).fieldIndex(), colIdx) = 1;

         SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);
         BoussinesqPerBetaCylGSystem::linearBlock(block, eqId, *fIt, nZ, kX_, kY_, Ra, Pr, Gamma, chi);
         tmp = Eigen::kroneckerProduct(blockMatrix, block.first);
         matrixRow.first += tmp;
         tmp = Eigen::kroneckerProduct(blockMatrix, block.second);
         matrixRow.second += tmp;

         colIdx++;
      }

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   DecoupledZSparse IBoussinesqPerBetaCylGScalarEquation::timeRow(FieldComponents::Spectral::Id comp, const int matIdx) const
   {
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Get mode indexes
      ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get radial wave number rescaled to box size
      MHDFloat kX_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1));
      // Get azimuthal wave number rescaled to box size
      MHDFloat kY_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0));

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = this->eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat Gamma = this->eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = this->eqParams().nd(NonDimensional::CHI);

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)), SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)),SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)));
      SparseMatrix  tmp(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx));

      // Create time row
      SparseMatrix   blockMatrix(this->couplingInfo(comp).nBlocks(),this->couplingInfo(comp).nBlocks());
      blockMatrix.insert(this->couplingInfo(comp).fieldIndex(), this->couplingInfo(comp).fieldIndex()) = 1;
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);
      BoussinesqPerBetaCylGSystem::timeBlock(block, eqId, nZ, kX_, kY_, Ra, Pr, Gamma, chi);
      tmp = Eigen::kroneckerProduct(blockMatrix, block.first);
      matrixRow.first += tmp;
      tmp = Eigen::kroneckerProduct(blockMatrix, block.second);
      matrixRow.second += tmp;

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   DecoupledZSparse IBoussinesqPerBetaCylGScalarEquation::boundaryRow(FieldComponents::Spectral::Id comp, const int matIdx) const
   {
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Get mode indexes
      ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get radial wave number rescaled to box size
      MHDFloat kX_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1));
      // Get azimuthal wave number rescaled to box size
      MHDFloat kY_ = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0));

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = this->eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat Gamma = this->eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = this->eqParams().nd(NonDimensional::CHI);

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)), SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)),SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)));
      SparseMatrix  tmp(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx));

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::field_iterator fIt;
      CouplingInformation::field_iterator_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix   blockMatrix(this->couplingInfo(comp).nBlocks(),this->couplingInfo(comp).nBlocks());
         blockMatrix.insert(this->couplingInfo(comp).fieldIndex(), colIdx) = 1;

         SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);
         BoussinesqPerBetaCylGSystem::boundaryBlock(block, eqId, *fIt, this->mspBcIds, nZ, kX_, kY_, Ra, Pr, Gamma, chi);
         tmp = Eigen::kroneckerProduct(blockMatrix, block.first);
         matrixRow.first += tmp;
         tmp = Eigen::kroneckerProduct(blockMatrix, block.second);
         matrixRow.second += tmp;

         colIdx++;
      }

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   void IBoussinesqPerBetaCylGScalarEquation::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      // Equation key
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Get X and Z dimensions
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Get the number of systems
      int nSystems = this->couplingInfo(FieldComponents::Spectral::SCALAR).nSystems();

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = this->eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat Gamma = this->eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = this->eqParams().nd(NonDimensional::CHI);

      // Boxscale
      MHDFloat boxYScale = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);
      MHDFloat boxXScale = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D);

      //
      // Initialise the quasi-inverse operators for the nonlinear terms
      //
      this->mNLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<SparseMatrix>()));
      std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::iterator qIt = this->mNLMatrices.find(FieldComponents::Spectral::SCALAR);
      qIt->second.reserve(nSystems);
      for(int i = 0; i < nSystems; ++i)
      {
         qIt->second.push_back(SparseMatrix());

         BoussinesqPerBetaCylGSystem::quasiInverse(qIt->second.back(), eqId, nZ);
      }

      //
      // Initialise the explicit linear operators
      //
      CouplingInformation::field_iterator fIt;
      CouplingInformation::field_iterator_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         std::vector<DecoupledZSparse> tmpMat;
         tmpMat.reserve(nSystems);

         bool isComplex = false;

         // Create matrices
         for(int i = 0; i < nSystems; ++i)
         {
            // Get mode indexes
            ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(i);

            // Get radial wave number rescaled to box size
            MHDFloat kX_ = boxXScale*static_cast<MHDFloat>(mode(1));
            // Get azimuthal wave number rescaled to box size
            MHDFloat kY_ = boxYScale*static_cast<MHDFloat>(mode(0));

            // Get linear block
            tmpMat.push_back(DecoupledZSparse());
            BoussinesqPerBetaCylGSystem::linearBlock(tmpMat.at(i), eqId, *fIt, nZ, kX_, kY_, Ra, Pr, Gamma, chi);

            // Explicit operator requires an additional minus sign
            tmpMat.at(i).first = -tmpMat.at(i).first;
            tmpMat.at(i).second = -tmpMat.at(i).second;

            isComplex = isComplex || (tmpMat.at(i).second.nonZeros() > 0);
         }

         // Create key
         std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(FieldComponents::Spectral::SCALAR, *fIt);

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
}
}
