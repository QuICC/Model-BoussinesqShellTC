/** \file IBoussinesqBetaCylGScalarEquation.cpp
 *  \brief Source of the implementation of the streamfunction equation in the 3DQG beta model
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/IBoussinesqBetaCylGScalarEquation.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGSystem.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IBoussinesqBetaCylGScalarEquation::IBoussinesqBetaCylGScalarEquation(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   IBoussinesqBetaCylGScalarEquation::~IBoussinesqBetaCylGScalarEquation()
   {
   }

   void IBoussinesqBetaCylGScalarEquation::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      // Equation key
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Get X and Z dimensions
      int nx = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nz = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get the number of systems
      int nSystems = this->couplingInfo(FieldComponents::Spectral::SCALAR).nSystems();

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = this->eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat Gamma = this->eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = this->eqParams().nd(NonDimensional::CHI);

      // Boxscale
      MHDFloat boxScale = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);

      //
      // Initialise the quasi-inverse operators for the nonlinear terms
      //
      this->mNLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<SparseMatrix>()));
      std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::iterator qIt = this->mNLMatrices.find(FieldComponents::Spectral::SCALAR);
      qIt->second.reserve(nSystems);
      for(int i = 0; i < nSystems; ++i)
      {
         qIt->second.push_back(SparseMatrix());

         BoussinesqBetaCylGSystem::quasiInverse(qIt->second.back(), eqId, nx, nz);
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
            MHDFloat k_ = boxScale*static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i));

            // Get linear block
            tmpMat.push_back(DecoupledZSparse());
            BoussinesqBetaCylGSystem::linearBlock(tmpMat.at(i), eqId, *fIt, nx, nz, k_, Ra, Pr, Gamma, chi);

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
