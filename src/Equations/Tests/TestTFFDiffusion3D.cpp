/** 
 * @file TestTFFDiffusion3D.cpp
 * @brief Source of the implementation of the TFF test equation for 3D diffusion
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
#include "Equations/Tests/TestTFFDiffusion3D.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "SpectralOperators/Tools/SpectralBoxTools.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestTFFDiffusion3D::TestTFFDiffusion3D(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   TestTFFDiffusion3D::~TestTFFDiffusion3D()
   {
   }

   void TestTFFDiffusion3D::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestTFFDiffusion3D::setCoupling()
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

   void TestTFFDiffusion3D::setRequirements()
   {
      // Add requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   DecoupledZSparse TestTFFDiffusion3D::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary) const
   {
      if(opId == IEquation::TIMEROW)
      { 
         return EigenSelector::timeRow(*this, compId, matIdx, hasBoundary);
      } else if(opId == IEquation::LINEARROW)
      {
         return EigenSelector::linearRow(*this, compId, matIdx, hasBoundary);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void linearBlock(const TestTFFDiffusion3D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
   {
      assert(eigs.size() == 2);
      // Rescale wave number to [-1, 1]
      MHDFloat k_ = eigs.at(0)/2.0;
      MHDFloat m_ = eigs.at(1)/2.0;

      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nX);

      EigenSelector::KZSum blocks;

      // Setup 3D diffusion
      if(fieldId.first == eq.name())
      {
         blocks.real() = Spectral::BoxTools::qLaplacian3D(spec1D, k_, m_, 2);

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }

      EigenSelector::constrainBlock(eq, compId, mat, fieldId, blocks, eigs, hasBoundary);
   }

   void timeBlock(const TestTFFDiffusion3D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
   {
      assert(eigs.size() == 2);

      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nX);

      EigenSelector::KZSum blocks;

      if(fieldId.first == eq.name())
      {
         blocks.real() = spec1D.qDiff(2,0);
      } else
      {
         throw Exception("Multiple field in time integration not implemented yet!");
      }

      EigenSelector::constrainBlock(eq, compId, mat, fieldId, blocks, eigs, hasBoundary);
   }

   void boundaryBlock(const TestTFFDiffusion3D& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx)
   {
      assert(eigs.size() == 2);

      if(fieldId.first == eq.name())
      {
         coeffs.push_back(1.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for boundary operator!");
      }
   }

}
}
