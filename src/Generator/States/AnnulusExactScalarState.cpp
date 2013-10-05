/** 
 * @file AnnulusExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in an annulus
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
#include "Generator/States/AnnulusExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   AnnulusExactScalarState::AnnulusExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(CONSTANT)
   {
   }

   AnnulusExactScalarState::~AnnulusExactScalarState()
   {
   }

   void AnnulusExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void AnnulusExactScalarState::setStateType(const AnnulusExactScalarState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void AnnulusExactScalarState::setCoupling()
   {
      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: first complex solver, complex solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, false, 0);

      // Set nonlinear flags: NO nonlinear term, NO quasi-inverse
      infoIt.first->second.setNonlinear(true, false);

      // Set source flags: has source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to itself
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Set mininal matrix coupling
      int nMat;
      ArrayI blockNs;
      ArrayI rhsCols;
      EigenSelector::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void AnnulusExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CONSTANT)
      {
         int nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, 1, 0.35);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array zGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               Array funcZ = 2.*(5.*(Math::PI/2)*(1+zGrid.array())).array().sin();
               Array funcR = 1 + rGrid.array() + rGrid.array().pow(3);
               MHDFloat funcTh = std::cos(thGrid(iTh));

               rNLComp.setProfile(funcZ*funcR(iR)*funcTh,iTh,iR);
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    MHDComplex AnnulusExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return MHDComplex(0,0);
    }

   void AnnulusExactScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   void AnnulusExactScalarState::createBoundaries(FieldComponents::Spectral::Id compId, const int matIdx)
   {
      EigenSelector::boundaryRow(*this, compId, matIdx);
   }

   DecoupledZSparse AnnulusExactScalarState::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::LINEARROW)
      {
         return EigenSelector::linearRow(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void AnnulusExactScalarState::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      Equations::quasiInverseBlock(*this, compId, mat);
   }

   void AnnulusExactScalarState::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
   {
      Equations::linearBlock(*this, compId, mat, fieldId, eigs);
   }

   void quasiInverseBlock(const AnnulusExactScalarState& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
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

   void linearBlock(const AnnulusExactScalarState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs)
   {
      assert(eigs.size() == 1);
      MHDFloat k = eigs.at(0);

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
      mat.rela() = Eigen::kroneckerProduct(spec3D.id(0), spec1D.id(0));

      // Prune matrices for safety
      mat.real().prune(1e-32);
      mat.imag().prune(1e-32);
   }

   void boundaryBlock(AnnulusExactScalarState& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs)
   {
      assert(eigs.size() == 1);
      MHDFloat k = eigs.at(0);

      std::vector<MHDFloat> coeffs;
      std::vector<Boundary::BCIndex>  bcIdx;

      coeffs.push_back(1.0);
      bcIdx.push_back(Boundary::BCIndex(EigenSelector::INDEPENDENT));

      coeffs.push_back(1.0);
      bcIdx.push_back(Boundary::BCIndex(EigenSelector::INDEPENDENT));

      EigenSelector::storeBoundaryCondition(eq, compId, fieldId, coeffs, bcIdx);
   }

}
}
