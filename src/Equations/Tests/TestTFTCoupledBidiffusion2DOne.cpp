/** 
 * @file TestTFTCoupledBidiffusion2DOne.cpp
 * @brief Source of the implementation of the TFT test equation for 2D diffusion (within 3D model)
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
#include "Equations/Tests/TestTFTCoupledBidiffusion2DOne.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "TypeSelectors/EquationToolsSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestTFTCoupledBidiffusion2DOne::TestTFTCoupledBidiffusion2DOne(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mSecondName(PhysicalNames::MAGNETIC)
   {
   }

   TestTFTCoupledBidiffusion2DOne::~TestTFTCoupledBidiffusion2DOne()
   {
   }

   void TestTFTCoupledBidiffusion2DOne::setIdentity(const PhysicalNames::Id first, const PhysicalNames::Id second)
   {
      // Set the name
      this->setName(first);

      // Set second name
      this->mSecondName = second;

      // Set the variable requirements
      this->setRequirements();
   }

   void TestTFTCoupledBidiffusion2DOne::setCoupling()
   {
      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: prognostic equation real solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, false, 0);

      // Set nonlinear flags: has nonlinear term, has quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to temperature equation (self)
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Equation is coupled to temperature equation (self)
      infoIt.first->second.addImplicitField(this->mSecondName, FieldComponents::Spectral::SCALAR);

      // Set mininal matrix coupling
      int nMat;
      ArrayI blockNs;
      ArrayI rhsCols;
      EquationToolsType::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void TestTFTCoupledBidiffusion2DOne::setRequirements()
   {
      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->mSecondName, FieldRequirement(true, true, false, false));
   }

   DecoupledZSparse TestTFTCoupledBidiffusion2DOne::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::TIMEROW)
      { 
         return EquationToolsType::timeRow(*this, compId, matIdx);
      } else if(opId == IEquation::LINEARROW)
      {
         return EquationToolsType::linearRow(*this, compId, matIdx);
      } else if(opId == IEquation::BOUNDARYROW)
      {
         return EquationToolsType::boundaryRow(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void linearBlock(const TestTFTCoupledBidiffusion2DOne& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs)
   {
      assert(eigs.size() == 1);
      MHDFloat k = eigs.at(0);

      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Initialise output matrices
      mat.real().resize(nX*nZ,nX*nZ);
      mat.imag().resize(nX*nZ,nX*nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Setup 3D diffusion
      if(fieldId.first == eq.name())
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         mat.real() = Eigen::kroneckerProduct(spec3D.id(0), Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 2));

      } else if(fieldId.first == eq.mSecondName)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         mat.real() = Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0));
         mat.real() = -mat.real();

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }

      // Prune matrices for safety
      mat.real().prune(1e-32);
      mat.imag().prune(1e-32);
   }

   void timeBlock(const TestTFTCoupledBidiffusion2DOne& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::vector<MHDFloat>& eigs)
   {
      assert(eigs.size() == 1);
      MHDFloat k = eigs.at(0);

      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Initialise output matrices
      mat.real().resize(nX*nZ,nX*nZ);
      mat.imag().resize(nX*nZ,nX*nZ);

      // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
      mat.real() = Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0));

      // Prune matrices for safety
      mat.real().prune(1e-32);
      mat.imag().prune(1e-32);
   }

   void boundaryBlock(const TestTFTCoupledBidiffusion2DOne& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs)
   {
      assert(eigs.size() == 1);
      MHDFloat k = eigs.at(0);

      int pX = 0;
      int pZ = 0;

      // Set boundary condition prefactors
      MHDFloat cX = 1.0;
      MHDFloat cZ = 0.0;

      // Compute boundary block operator
      EquationToolsType::boundaryBlock(eq, FieldComponents::Spectral::SCALAR, mat, fieldId, pX, pZ, cX, cZ);
   }

}
}
