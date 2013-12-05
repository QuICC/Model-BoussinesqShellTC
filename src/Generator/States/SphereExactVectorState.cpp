/** 
 * @file SphereExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact vector solution in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <tr1/cmath>

// External includes
//

// Class include
//
#include "Generator/States/SphereExactVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   SphereExactVectorState::SphereExactVectorState(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mTypeId(CONSTANT)
   {
   }

   SphereExactVectorState::~SphereExactVectorState()
   {
   }

   void SphereExactVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set toroidal and poloidal components
      this->mSpectralIds.push_back(FieldComponents::Spectral::ONE);
      this->mSpectralIds.push_back(FieldComponents::Spectral::TWO);

      // Set the variable requirements
      this->setRequirements();
   }

   void SphereExactVectorState::setStateType(const SphereExactVectorState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void SphereExactVectorState::setHarmonicOptions(const std::vector<std::tr1::tuple<int,int,MHDComplex> >& rModes, const std::vector<std::tr1::tuple<int,int,MHDComplex> >& tModes, const std::vector<std::tr1::tuple<int,int,MHDComplex> >& pModes)
   {
      this->mRSHModes = rModes;
      this->mTSHModes = tModes;
      this->mPSHModes = pModes;
   }

   void SphereExactVectorState::setCoupling()
   {
      SpectralComponent_range specRange = this->spectralRange();
      for(SpectralComponent_iterator it = specRange.first; it != specRange.second; ++it)
      {
         // Initialise coupling information
         std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
         infoIt = this->mCouplingInfos.insert(std::make_pair(*it,CouplingInformation()));
         SpectralFieldId eqId = std::make_pair(this->name(), *it);

         // General setup: first complex solver, complex solver, start from m = 0
         infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, false, 0);

         // Set nonlinear flags: NO nonlinear term, NO quasi-inverse
         infoIt.first->second.setNonlinear(true, false);

         // Set source flags: has source term
         infoIt.first->second.setSource(false);

         // Set index type: use MODE
         infoIt.first->second.setIndexType(CouplingInformation::MODE);

         // Equation is coupled to itself
         infoIt.first->second.addImplicitField(eqId.first, *it);

         // Set mininal matrix coupling
         int nMat = 0;
         ArrayI blockNs;
         ArrayI rhsCols;
         EigenSelector::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
         infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

         // Sort implicit fields
         infoIt.first->second.sortImplicitFields(eqId.first, *it);
      }
   }

   void SphereExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == CONSTANT)
      {
         rNLComp.rData().setConstant(42);
      } else if(this->mTypeId == HARMONIC)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array funcPh(nPh);
         MHDFloat funcR(nR);
         MHDFloat funcTh;
         typedef std::vector<std::tr1::tuple<int,int,MHDComplex> >::const_iterator ModeIt;
         std::pair<ModeIt, ModeIt>  modeRange;
         if(compId == FieldComponents::Physical::ONE)
         {
            modeRange.first = this->mRSHModes.begin();
            modeRange.second = this->mRSHModes.end();
         } else if(compId == FieldComponents::Physical::TWO)
         {
            modeRange.first = this->mTSHModes.begin();
            modeRange.second = this->mTSHModes.end();
         } else
         {
            modeRange.first = this->mPSHModes.begin();
            modeRange.second = this->mPSHModes.end();
         }
         rNLComp.rData().setConstant(0);
         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               for(ModeIt it = modeRange.first; it != modeRange.second; ++it)
               {
                  int l = std::tr1::get<0>(*it);
                  int m = std::tr1::get<1>(*it);
                  MHDFloat re = std::tr1::get<2>(*it).real();
                  MHDFloat im = std::tr1::get<2>(*it).imag();

                  funcR = std::pow(rGrid(iR),l);

                  // Spherical harmonic Y_l^m
                  funcPh = re*(static_cast<MHDFloat>(m)*phGrid).array().cos() + im*(static_cast<MHDFloat>(m)*phGrid).array().sin();
                  funcTh = std::tr1::sph_legendre(l,m, thGrid(iTh));
                  rNLComp.addProfile(funcPh*funcR*funcTh,iTh,iR);
               }
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    MHDComplex SphereExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      return MHDComplex(0,0);
    }

   void SphereExactVectorState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, false, false));
   }

   DecoupledZSparse SphereExactVectorState::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary) const
   {
      if(opId == IEquation::LINEARROW)
      {
         return EigenSelector::linearRow(*this, compId, matIdx, hasBoundary);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void SphereExactVectorState::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      Equations::quasiInverseBlock(*this, compId, mat);
   }

   void SphereExactVectorState::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
   {
      Equations::linearBlock(*this, compId, mat, fieldId, eigs, false);
   }

   void quasiInverseBlock(const SphereExactVectorState& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nR);

      // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B) => out = A(i,j)*B)
      mat = spec1D.id(0);

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const SphereExactVectorState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
   {
      assert(eigs.size() == 2);
      MHDFloat l = eigs.at(0);
      MHDFloat m = eigs.at(1);

      // Get X and Z dimensions
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nR);

      // Initialise output matrices
      mat.real().resize(nR,nR);
      mat.imag().resize(nR,nR);

      // Build linear operator (kronecker(A,B) => out = A(i,j)*B)
      mat.real() = spec1D.id(0);

      // Prune matrices for safety
      mat.real().prune(1e-32);
      mat.imag().prune(1e-32);
   }

   void boundaryBlock(SphereExactVectorState& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx)
   {
      assert(eigs.size() == 2);
      MHDFloat l = eigs.at(0);
      MHDFloat m = eigs.at(1);

      coeffs.push_back(1.0);
      bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));
   }

}
}