/** 
 * @file ExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution
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
#include "Generator/States/ExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   ExactScalarState::ExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(CONSTANT), mSineA(2), mSineN(2)
   {
   }

   ExactScalarState::~ExactScalarState()
   {
   }

   void ExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ExactScalarState::setStateType(const ExactScalarState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void ExactScalarState::setSineOptions(const MHDFloat aX, const MHDFloat kX, const MHDFloat aZ, const MHDFloat kZ)
   {
      this->mSineA(0) = aX;
      this->mSineA(1) = aZ;

      this->mSineN(0) = kX;
      this->mSineN(1) = kZ;
   }

   void ExactScalarState::setCoupling()
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
      int nMat = 0;
      ArrayI blockNs;
      ArrayI rhsCols;
      EigenSelector::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void ExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CONSTANT)
      {
         rNLComp.rData().setConstant(42.0);
      } else if(this->mTypeId == SINESINE)
      {
         int nX = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int nY = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array xGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nX);
         Array yGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nY);
         Array zGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

         for(int iX = 0; iX < nX; ++iX)
         {
            for(int iY = 0; iY < nY; ++iY)
            {
               Array sinZ = this->mSineA(1)*(this->mSineN(1)*(Math::PI/2)*(1+zGrid.array())).array().sin();
               Array sinX = this->mSineA(0)*(this->mSineN(0)*(Math::PI/2)*(1+xGrid.array())).array().sin();
               MHDFloat yVal = std::cos(yGrid(iY));
               //MHDFloat yVal = 1.0;

               rNLComp.setProfile(sinZ*sinX(iX)*yVal,iY,iX);
            }
         }
      } else if(this->mTypeId == SINECOSINE)
      {
         int nX = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int nY = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array xGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nX);
         Array yGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nY);
         Array zGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

         for(int iX = 0; iX < nX; ++iX)
         {
            for(int iY = 0; iY < nY; ++iY)
            {
               Array sinZ = this->mSineA(1)*(this->mSineN(1)*(Math::PI/2)*(1+zGrid.array())).array().cos();
               Array sinX = this->mSineA(0)*(this->mSineN(0)*(Math::PI/2)*(1+xGrid.array())).array().sin();
               MHDFloat yVal = std::cos(yGrid(iY));
               //MHDFloat yVal = 1.0;

               rNLComp.setProfile(sinZ*sinX(iX)*yVal,iY,iX);
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    MHDComplex ExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return MHDComplex(0,0);
    }

   void ExactScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

   DecoupledZSparse ExactScalarState::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary) const
   {
      if(opId == IEquation::LINEARROW)
      {
         return EigenSelector::linearRow(*this, compId, matIdx, hasBoundary);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void ExactScalarState::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      Equations::quasiInverseBlock(*this, compId, mat);
   }

   void ExactScalarState::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
   {
      Equations::linearBlock(*this, compId, mat, fieldId, eigs, false);
   }

   void quasiInverseBlock(const ExactScalarState& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nX);
      Spectral::OperatorSelector<Dimensions::Simulation::SIM3D>::Type spec3D(nZ);

      EigenSelector::KRSum blocks;
      EigenSelector::KRProduct kProduct(SparseMatrix(nX*nZ,nX*nZ),SparseMatrix(nX*nZ,nX*nZ));

      std::tr1::get<0>(kProduct) = spec1D.id(0);
      std::tr1::get<1>(kProduct) = spec3D.id(0);
      blocks.push_back(kProduct);

      EigenSelector::computeKSum(mat, blocks);
   }

   void linearBlock(const ExactScalarState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
   {
      assert(eigs.size() == 1);
      MHDFloat k = eigs.at(0);

      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nX);
      Spectral::OperatorSelector<Dimensions::Simulation::SIM3D>::Type spec3D(nZ);

      EigenSelector::KZSum blocks;
      EigenSelector::KZProduct kProduct(DecoupledZSparse(nX*nZ,nX*nZ),DecoupledZSparse(nX*nZ,nX*nZ));

      std::tr1::get<0>(kProduct).real() = spec1D.id(0);
      std::tr1::get<1>(kProduct).real() = spec3D.id(0);
      blocks.push_back(kProduct);

      EigenSelector::constrainBlock(eq, compId, mat, fieldId, blocks, eigs, hasBoundary);
   }

   void boundaryBlock(const ExactScalarState& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx)
   {  
      assert(eigs.size() == 1);
      MHDFloat k = eigs.at(0);

      coeffs.push_back(1.0);
      bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

      coeffs.push_back(1.0);
      bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));
   }

}
}
