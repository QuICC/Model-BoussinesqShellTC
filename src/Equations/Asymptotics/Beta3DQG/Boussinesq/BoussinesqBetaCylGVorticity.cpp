/** \file BoussinesqBetaCylGVorticity.cpp
 *  \brief Source of the implementation of the vorticity computation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGVorticity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Tools/Equation1DEigenTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBetaCylGVorticity::BoussinesqBetaCylGVorticity(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBetaCylGVorticity::~BoussinesqBetaCylGVorticity()
   {
   }

   void BoussinesqBetaCylGVorticity::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      this->initSpectralMatrices1DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
   }

   void BoussinesqBetaCylGVorticity::setCoupling()
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

      // General setup: trivial equation, complex solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, false, 0);

      // set nonlinear flags: NO nonlinear term, NO quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to vorticity equation (self)
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Equation has explicit temperature
      infoIt.first->second.addExplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR);

      // Set sizes of blocks and matrices
      ArrayI blockNs(nY);
      blockNs.setConstant(nX*nZ);
      ArrayI rhsCols(nY);
      rhsCols.setConstant(1);
      infoIt.first->second.setSizes(nY, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void BoussinesqBetaCylGVorticity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Doesn't require a nonlinear computation
      assert(false);
   }

   void BoussinesqBetaCylGVorticity::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VORTICITY);

      // Set vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITY, FieldRequirement(true, true, true, true));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, false));
   }

   DecoupledZSparse BoussinesqBetaCylGVorticity::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::LINEARROW)
      {
         return Equation1DEigenTools::linearRow(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void BoussinesqBetaCylGVorticity::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      // Quasi inverse is not required without nonlinear computations
      assert(false);
   }

   void BoussinesqBetaCylGVorticity::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      linearBlock(*this, compId, mat, fieldId, k);
   }

   void linearBlock(const BoussinesqBetaCylGVorticity& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
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

      /// - Streamfunction: \f$ \zeta = \nabla^2 \psi\f$
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         mat.first = Eigen::kroneckerProduct(spec3D.id(0), Spectral::PeriodicOperator::laplacian2D(spec1D, k_, 0));

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

}
}