/** \file TestTFTDiffusion2D.cpp
 *  \brief Source of the implementation of the TFT test equation for 2D diffusion (within 3D model)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Tests/TestTFTDiffusion2D.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Tools/Equation1DEigenTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestTFTDiffusion2D::TestTFTDiffusion2D(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   TestTFTDiffusion2D::~TestTFTDiffusion2D()
   {
   }

   void TestTFTDiffusion2D::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestTFTDiffusion2D::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      this->initSpectralMatrices1DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
   }

   void TestTFTDiffusion2D::setCoupling()
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

      // General setup: prognostic equation real solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, false, 0);

      // Set nonlinear flags: has nonlinear term, has quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to temperature equation (self)
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

   void TestTFTDiffusion2D::setRequirements()
   {
      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   DecoupledZSparse TestTFTDiffusion2D::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::TIMEROW)
      { 
         return Equation1DEigenTools::timeRow(*this, compId, matIdx);
      } else if(opId == IEquation::LINEARROW)
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

   void linearBlock(const TestTFTDiffusion2D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
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

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Setup 3D diffusion
      if(fieldId.first == eq.name())
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         mat.first = Eigen::kroneckerProduct(spec3D.id(0), Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 2));

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void timeBlock(const TestTFTDiffusion2D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat k)
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
      mat.first = Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0));

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void boundaryBlock(const TestTFTDiffusion2D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      int pX = 0;
      int pZ = 0;

      // Set boundary condition prefactors
      MHDFloat cX = 1.0;
      MHDFloat cZ = 0.0;

      // Compute boundary block operator
      Equation1DEigenTools::boundaryBlock1DEigen(eq, FieldComponents::Spectral::SCALAR, mat, fieldId, pX, pZ, cX, cZ);
   }

}
}
