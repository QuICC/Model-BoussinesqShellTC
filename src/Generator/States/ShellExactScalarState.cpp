/** \file ShellExactScalarState.cpp
 *  \brief Source of the implementation of the equation to generate an exact scalar solution in a spherical shell
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Generator/States/ShellExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"

#include <iostream>
namespace GeoMHDiSCC {

namespace Equations {

   ShellExactScalarState::ShellExactScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mTypeId(CONSTANT)
   {
   }

   ShellExactScalarState::~ShellExactScalarState()
   {
   }

   void ShellExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ShellExactScalarState::setStateType(const ShellExactScalarState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void ShellExactScalarState::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      this->initSpectralMatrices1DPeriodic(spBcIds);
   }

   void ShellExactScalarState::setCoupling()
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

   void ShellExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CONSTANT)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, 1, 0.35);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array funcPh(nPh);
         MHDFloat funcR(nR);
         MHDFloat funcTh;
         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               funcR = 1.0;

               // Spherical harmonic Y_0^0
               funcPh.setConstant(1);
               funcTh = 0.5*std::sqrt(1./MathConstants::PI);
               rNLComp.setProfile(funcPh*funcR*funcTh,iTh,iR);

               // Spherical harmonic Y_1^0
               funcPh.setConstant(1);
               funcTh = 0.5*std::sqrt(3./(MathConstants::PI))*std::cos(thGrid(iTh));
               rNLComp.addProfile(funcPh*funcR*funcTh,iTh,iR);

               // Spherical harmonic Y_1^1
               funcPh = phGrid.array().sin();
               funcTh = -0.5*std::sqrt(3./(2.*MathConstants::PI))*std::sin(thGrid(iTh));
               rNLComp.addProfile(funcPh*funcR*funcTh,iTh,iR);

               // Spherical harmonic Y_2^0
               funcPh.setConstant(1);
               funcTh = 0.25*std::sqrt(5./MathConstants::PI)*(3.*std::pow(std::cos(thGrid(iTh)),2)-1);
               rNLComp.addProfile(funcPh*funcR*funcTh,iTh,iR);

               // Spherical harmonic Y_2^1
               funcPh = phGrid.array().cos();
               funcTh = -0.5*std::sqrt(15./(2.*MathConstants::PI))*std::sin(thGrid(iTh))*std::cos(thGrid(iTh));
               rNLComp.addProfile(funcPh*funcR*funcTh,iTh,iR);

               // Spherical harmonic Y_2^2
               funcPh = (2.*phGrid).array().cos();
               funcTh = 0.25*std::sqrt(15./(2.*MathConstants::PI))*std::pow(std::sin(thGrid(iTh)),2);
               rNLComp.addProfile(funcPh*funcR*funcTh,iTh,iR);
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    MHDComplex ShellExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return MHDComplex(0,0);
    }

   void ShellExactScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   DecoupledZSparse ShellExactScalarState::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::LINEARROW)
      {
         throw Exception("Operators for 2D eigen directions not implemented yet");
      } else if(opId == IEquation::BOUNDARYROW)
      {
         throw Exception("Operators for 2D eigen directions not implemented yet");
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void ShellExactScalarState::setQuasiInverse(SparseMatrix& mat) const
   {
      Equations::quasiInverseBlock(*this, mat);
   }

   void ShellExactScalarState::setExplicitLinearBlock(DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      Equations::linearBlock(*this, mat, fieldId, k);
   }

   void quasiInverseBlock(const ShellExactScalarState& eq, SparseMatrix& mat)
   {
      throw Exception("Operators for 2D eigen directions not implemented yet");
   }

   void linearBlock(const ShellExactScalarState& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      throw Exception("Operators for 2D eigen directions not implemented yet");
   }

   void boundaryBlock(const ShellExactScalarState& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      // Compute boundary block operator
      throw Exception("Operators for 2D eigen directions not implemented yet");
   }

}
}
