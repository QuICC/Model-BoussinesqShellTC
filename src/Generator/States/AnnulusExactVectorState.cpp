/** 
 * @file AnnulusExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact vector solution in an annulus
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
#include "Generator/States/AnnulusExactVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   AnnulusExactVectorState::AnnulusExactVectorState(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mTypeId(CONSTANT)
   {
   }

   AnnulusExactVectorState::~AnnulusExactVectorState()
   {
   }

   void AnnulusExactVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      this->mSpectralIds.push_back(FieldComponents::Spectral::ONE);
      this->mSpectralIds.push_back(FieldComponents::Spectral::TWO);

      // Set the variable requirements
      this->setRequirements();
   }

   void AnnulusExactVectorState::setStateType(const AnnulusExactVectorState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void AnnulusExactVectorState::setCoupling()
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

   void AnnulusExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == CONSTANT)
      {
         int nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, 1.0, 0.35);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array zGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               for(int iZ = 0; iZ < nZ; ++iZ)
               {
                  rNLComp.setPoint(1.0, iZ,iTh,iR);
               }
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    MHDComplex AnnulusExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      return MHDComplex(0,0);
    }

   void AnnulusExactVectorState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, false, false));
   }

   DecoupledZSparse AnnulusExactVectorState::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary) const
   {
      if(opId == IEquation::LINEARROW)
      {
         return EigenSelector::linearRow(*this, compId, matIdx, hasBoundary);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void AnnulusExactVectorState::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      Equations::quasiInverseBlock(*this, compId, mat);
   }

   void AnnulusExactVectorState::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
   {
      Equations::linearBlock(*this, compId, mat, fieldId, eigs, false);
   }

   void quasiInverseBlock(const AnnulusExactVectorState& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nX);
      Spectral::OperatorSelector<Dimensions::Simulation::SIM3D>::Type spec3D(nZ);

      // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B) => out = A(i,j)*B)
      mat = Eigen::kroneckerProduct(spec3D.id(0), spec1D.id(0));

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const AnnulusExactVectorState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
   {
      assert(eigs.size() == 1);
      MHDFloat m = eigs.at(0);

      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nX);
      Spectral::OperatorSelector<Dimensions::Simulation::SIM3D>::Type spec3D(nZ);

      // Initialise output matrices
      mat.real().resize(nX*nZ,nX*nZ);
      mat.imag().resize(nX*nZ,nX*nZ);

      // Build linear operator (kronecker(A,B) => out = A(i,j)*B)
      mat.real() = Eigen::kroneckerProduct(spec3D.id(0), spec1D.id(0));

      // Prune matrices for safety
      mat.real().prune(1e-32);
      mat.imag().prune(1e-32);
   }

   void boundaryBlock(AnnulusExactVectorState& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx)
   {
      assert(eigs.size() == 1);
      MHDFloat m = eigs.at(0);

      coeffs.push_back(1.0);
      bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

      coeffs.push_back(1.0);
      bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));
   }

}
}
