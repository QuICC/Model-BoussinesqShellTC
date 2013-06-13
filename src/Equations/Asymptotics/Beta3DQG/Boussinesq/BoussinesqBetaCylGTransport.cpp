/** \file BoussinesqBetaCylGTransport.cpp
 *  \brief Source of the implementation of the transport equation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGSystem.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBetaCylGTransport::BoussinesqBetaCylGTransport(SharedEquationParameters spEqParams)
      : IBoussinesqBetaCylGScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBetaCylGTransport::~BoussinesqBetaCylGTransport()
   {
   }

   void BoussinesqBetaCylGTransport::setCoupling()
   {
      // Get X dimension
      int nX = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get Y dimension
      int nY = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Initialise coupling information
      this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: first real solver, real solver, start from m = 0
      this->mCouplingInfo.setGeneral(0, false, 0);

      // 
      //  WARNING: the order is important
      //

      // Equation is coupled to temperature equation (self)
      this->mCouplingInfo.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR, true);

      // Equation has explicit streamfunction
      this->mCouplingInfo.addExplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR);

      // Set sizes of blocks and matrices
      ArrayI blockNs(nY);
      blockNs.setConstant(nX*nZ);
      ArrayI rhsCols(nY);
      rhsCols.setConstant(1);
      this->mCouplingInfo.setSizes(nY, blockNs, rhsCols); 
   }

   void BoussinesqBetaCylGTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T}\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 0.0);
   }

   void BoussinesqBetaCylGTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, true));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

   DecoupledZSparse linearRow(const BoussinesqBetaCylGTransport& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      return linearRow1DPeriodic(eq, compId, matIdx);
   }

   DecoupledZSparse timeRow(const BoussinesqBetaCylGTransport& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      return timeRow1DPeriodic(eq, compId, matIdx);
   }

   DecoupledZSparse boundaryRow(const BoussinesqBetaCylGTransport& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      return boundaryRow1DPeriodic(eq, compId, matIdx);
   }

   void linearBlock(const BoussinesqBetaCylGTransport& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int nX, const int nZ, const MHDFloat k)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Initialise output matrices
      mat.first.resize(nX*nZ,nX*nZ);
      mat.second.resize(nX*nZ,nX*nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = eq.eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = eq.eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat Gamma = eq.eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = eq.eqParams().nd(NonDimensional::CHI);

      /// - Streamfunction : \f$ i \frac{k}{2} \left(D_x^{-2} \otimes I_Z\right) \f$
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         SparseMatrix tmp = k_*spec3D.id(0);
         Eigen::kroneckerProduct(tmp, spec1D.qDiff(2,0), mat.second);

         /// - Vertical velocity : \f$ \left(0_x \otimes 0_Z\right) \f$
      } else if(fieldId.first == PhysicalNames::VELOCITYZ)
      {
         //
         // Nothing to do for an empty sparse matrix
         //

         /// - Temperature : \f$ \frac{1}{Pr}\left(D_x^{-2}\nabla_\perp^{2} \otimes I_Z\right) \f$
      } else if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         SparseMatrix tmp = (1./Pr)*spec3D.id(0);
         Eigen::kroneckerProduct(tmp, Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 2), mat.first);

         // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void timeBlock(const BoussinesqBetaCylGTransport& eq, DecoupledZSparse& mat, const int nX, const int nZ, const MHDFloat k)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Initialise output matrices
      mat.first.resize(nX*nZ,nX*nZ);
      mat.second.resize(nX*nZ,nX*nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
      Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0), mat.first);

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void boundaryBlock(const BoussinesqBetaCylGTransport& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const SharedSimulationBoundary spBcIds, const int nX, const int nZ, const MHDFloat k)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Create spectral boundary operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::BcType bound1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::BcType bound3D(nZ);

      // Initialise output matrices
      mat.first.resize(nX*nZ,nX*nZ);
      mat.second.resize(nX*nZ,nX*nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Storage for the boundary quasi-inverses
      SparseMatrix q1D;
      SparseMatrix q3D;

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

}
}
