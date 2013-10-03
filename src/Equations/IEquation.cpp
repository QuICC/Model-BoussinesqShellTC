/** 
 * @file IEquation.cpp
 * @brief Source of building block for the implementation of a time dependend evolution equation
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
#include "Equations/IEquation.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/EquationToolsSelector.hpp"

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

   void IEquation::initSpectralMatricesComponent(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId)
   {
      // Get the number of systems
      int nSystems = this->couplingInfo(compId).nSystems();

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
            std::vector<MHDFloat> eigs = EquationToolsType::getEigs(*this, i);

            // Get linear block
            tmpMat.push_back(DecoupledZSparse());
            this->setExplicitLinearBlock(compId, tmpMat.at(i), *fIt, eigs);

            // Explicit operator requires an additional minus sign
            tmpMat.at(i).real() = -tmpMat.at(i).real();
            tmpMat.at(i).imag() = -tmpMat.at(i).imag();

            isComplex = isComplex || (tmpMat.at(i).imag().nonZeros() > 0);
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
               SparseMatrixZ tmp = tmpMat.at(i).real().cast<MHDComplex>() + MathConstants::cI*tmpMat.at(i).imag();
               this->mLZMatrices.find(key)->second.push_back(tmp);
            }
         } else
         {
            this->mLDMatrices.insert(std::make_pair(key, std::vector<SparseMatrix>()));
            this->mLDMatrices.find(key)->second.back().reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               this->mLDMatrices.find(key)->second.push_back(SparseMatrix());

               this->mLDMatrices.find(key)->second.back().swap(tmpMat.at(i).real());
            }
         }
      }
   }

   void IEquation::setQuasiInverse(FieldComponents::Spectral::Id comp, SparseMatrix &mat) const
   {
      throw Exception("setQuasiInverse: dummy implementation was called!");
   }

   void IEquation::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
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

   void quasiInverseBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      throw Exception("quasiInverseBlock: dummy implementation called!");
   }

   void timeBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::vector<MHDFloat>& k)
   {
      throw Exception("timeBlock: dummy implementation called!");
   }

   void linearBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& k)
   {
      throw Exception("linearBlock: dummy implementation called!");
   }

   void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& k)
   {
      throw Exception("boundaryBlock: dummy implementation called!");
   }
}
}
