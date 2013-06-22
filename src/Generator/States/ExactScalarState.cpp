/** \file ExactScalarState.cpp
 *  \brief Source of the implementation of the equation to generate an exact scalar solution
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
#include "TypeSelectors/SpectralSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "Equations/Tools/Equation1DEigenTools.hpp"

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

   void ExactScalarState::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      this->initSpectralMatrices1DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
   }

   void ExactScalarState::setCoupling()
   {
      // Get X dimension
      int nX = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get Y dimension
      int nY = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

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

      // Set sizes of blocks and matrices
      ArrayI blockNs(nY);
      blockNs.setConstant(nX*nZ);
      ArrayI rhsCols(nY);
      rhsCols.setConstant(1);
      infoIt.first->second.setSizes(nY, blockNs, rhsCols); 

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
      } else if(this->mTypeId == SINE)
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
               Array sinZ = this->mSineA(1)*(this->mSineN(1)*(MathConstants::PI/2)*(1+zGrid.array())).array().sin();
               Array sinX = this->mSineA(0)*(this->mSineN(0)*(MathConstants::PI/2)*(1+xGrid.array())).array().sin();
               //MHDFloat yVal = std::cos(yGrid(iY));
               MHDFloat yVal = 1.0;

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
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   DecoupledZSparse ExactScalarState::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::LINEARROW)
      {
         return Equation1DEigenTools::linearRow(*this, compId, matIdx);
      } else if(opId == IEquation::BOUNDARYROW)
      {
         return Equation1DEigenTools::boundaryRow(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void ExactScalarState::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      Equations::quasiInverseBlock(*this, compId, mat);
   }

   void ExactScalarState::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      Equations::linearBlock(*this, compId, mat, fieldId, k);
   }

   void quasiInverseBlock(const ExactScalarState& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
      Eigen::kroneckerProduct(spec3D.id(0), spec1D.id(0), mat);

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const ExactScalarState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Initialise output matrices
      mat.first.resize(nX*nZ,nX*nZ);
      mat.second.resize(nX*nZ,nX*nZ);

      // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
      Eigen::kroneckerProduct(spec3D.id(0), spec1D.id(0), mat.first);

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void boundaryBlock(const ExactScalarState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      int pX = 0;
      int pZ = 0;

      // Set boundary condition prefactors
      MHDFloat cX = 1.0;
      MHDFloat cZ = 1.0;

      // Compute boundary block operator
      Equation1DEigenTools::boundaryBlock1DEigen(eq, FieldComponents::Spectral::SCALAR, mat, fieldId, pX, pZ, cX, cZ);
   }

}
}
