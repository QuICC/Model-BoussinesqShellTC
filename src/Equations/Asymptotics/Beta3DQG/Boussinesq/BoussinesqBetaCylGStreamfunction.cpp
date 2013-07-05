/** \file BoussinesqBetaCylGStreamfunction.cpp
 *  \brief Source of the implementation of the streamfunction equation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGStreamfunction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Tools/Equation1DEigenTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBetaCylGStreamfunction::BoussinesqBetaCylGStreamfunction(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBetaCylGStreamfunction::~BoussinesqBetaCylGStreamfunction()
   {
   }

   void BoussinesqBetaCylGStreamfunction::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      this->initSpectralMatrices1DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
   }

   void BoussinesqBetaCylGStreamfunction::setCoupling()
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

      // General setup: prognostic equation, complex solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, true, 0);

      // set nonlinear flags: has nonlinear term, has quasi-inverse
      infoIt.first->second.setNonlinear(true, true);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to streamfunction equation (self)
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Equation is coupled to vertical velocity equation
      infoIt.first->second.addImplicitField(PhysicalNames::VELOCITYZ,FieldComponents::Spectral::SCALAR);

      infoIt.first->second.addImplicitField(PhysicalNames::TEMPERATURE,FieldComponents::Spectral::SCALAR);

      // Equation has explicit temperature
      //infoIt.first->second.addExplicitField(PhysicalNames::TEMPERATURE,FieldComponents::Spectral::SCALAR);

      // Set sizes of blocks and matrices
      ArrayI blockNs(nY);
      blockNs.setConstant(nX*nZ);
      ArrayI rhsCols(nY);
      rhsCols.setConstant(1);
      infoIt.first->second.setSizes(nY, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void BoussinesqBetaCylGStreamfunction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\nabla^2_{\perp}\psi\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->unknown().dom(0).grad(), this->scalar(PhysicalNames::VORTICITYZ).dom(0).grad(), 1.0);
   }

   void BoussinesqBetaCylGStreamfunction::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::STREAMFUNCTION);

      // Set streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, true));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, false, false));

      // Add temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, false, false, false));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, false, true));
   }

   DecoupledZSparse BoussinesqBetaCylGStreamfunction::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
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

   void BoussinesqBetaCylGStreamfunction::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      quasiInverseBlock(*this, compId, mat);
   }

   void BoussinesqBetaCylGStreamfunction::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      linearBlock(*this, compId, mat, fieldId, k);
   }

   void quasiInverseBlock(const BoussinesqBetaCylGStreamfunction& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      /// - Streamfunction equation: \f$ \left(D_x^{-4} \otimes D_Z^{-1}\right) \f$
      // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
      Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(4,0), mat);

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const BoussinesqBetaCylGStreamfunction& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
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

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = eq.eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = eq.eqParams().nd(NonDimensional::PRANDTL);

      /// - Streamfunction : \f$ \left(D_x^{-4}\nabla_\perp^{4} \otimes D_Z^{-1}\right) \f$
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), Spectral::PeriodicOperator::qBilaplacian2D(spec1D, k_, 4), mat.first);

         /// - Vertical velocity : \f$ \left(D_x^{-4} \otimes I_z^{-1}\right) \f$
      } else if(fieldId.first == PhysicalNames::VELOCITYZ)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.id(1), spec1D.qDiff(4,0), mat.first);

         /// - Temperature : \f$ i \frac{k}{2}\frac{1}{16}\frac{Ra}{Pr}\left( D_x^{-4} \otimes D_Z^{-1}\right) \f$
      } else if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         SparseMatrix tmp = k_*(1.0/16.)*(Ra/Pr)*spec3D.qDiff(1,0);
         Eigen::kroneckerProduct(tmp, spec1D.qDiff(4,0), mat.second);

         // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void timeBlock(const BoussinesqBetaCylGStreamfunction& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat k)
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

      // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
      Eigen::kroneckerProduct(spec3D.qDiff(1,0), Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 4), mat.first);

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void boundaryBlock(const BoussinesqBetaCylGStreamfunction& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Get Physical parameters Gamma, chi
      MHDFloat Gamma = eq.eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = eq.eqParams().nd(NonDimensional::CHI);

      int pX = 1;
      int pZ = 0;

      MHDFloat cX;
      MHDFloat cZ;

      // Boundary condition for the streamfunction
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Set boundary condition prefactors
         cX = 1.0;
         cZ = k_*std::tan((MathConstants::PI/180.)*chi)/Gamma;

      // Boundary condition for the vertical velocity
      } else if(fieldId.first == PhysicalNames::VELOCITYZ)
      {
         // Set boundary condition prefactors
         cX = 1.0;
         cZ = 1.0;

      // Boundary condition for the temperature
      } else if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         // Set boundary condition prefactors
         cX = 0.0;
         cZ = 0.0;

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for boundary operator!");
      }

      // Compute boundary block operator
      Equation1DEigenTools::boundaryBlock1DEigen(eq, FieldComponents::Spectral::SCALAR, mat, fieldId, pX, pZ, cX, cZ);
   }

}
}
