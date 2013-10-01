/** 
 * @file BoussinesqBetaCylGTransport.cpp
 * @brief Source of the implementation of the transport equation in the 3DQG beta model
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "TypeSelectors/EquationToolsSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBetaCylGTransport::BoussinesqBetaCylGTransport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBetaCylGTransport::~BoussinesqBetaCylGTransport()
   {
   }

   void BoussinesqBetaCylGTransport::setCoupling()
   {
      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: prognostic equation, real solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, true, 0);

      // Set nonlinear flags: has nonlinear term, has quasi-inverse
      infoIt.first->second.setNonlinear(true, true);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to temperature equation (self)
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      infoIt.first->second.addImplicitField(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);

      infoIt.first->second.addImplicitField(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);

      // Equation has explicit streamfunction
      //infoIt.first->second.addExplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR);

      // Set mininal matrix coupling
      int nMat;
      ArrayI blockNs;
      ArrayI rhsCols;
      EquationToolsType::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void BoussinesqBetaCylGTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T}\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);
   }

   void BoussinesqBetaCylGTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, true));

      // Add streamfunction to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

   DecoupledZSparse BoussinesqBetaCylGTransport::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
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

   void BoussinesqBetaCylGTransport::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      quasiInverseBlock(*this, compId, mat);
   }

   void BoussinesqBetaCylGTransport::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      linearBlock(*this, compId, mat, fieldId, k);
   }

   void quasiInverseBlock(const BoussinesqBetaCylGTransport& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      /// - Transport equation: \f$ \left( D_x^{-2} \otimes I_Z\right) \f$
      // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
      mat = Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0));

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const BoussinesqBetaCylGTransport& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
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

      // Get Physical parameters Pr
      MHDFloat Pr = eq.eqParams().nd(NonDimensional::PRANDTL);

      /// - Streamfunction : \f$ i \frac{k}{2} \left(D_x^{-2} \otimes I_Z\right) \f$
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         SparseMatrix tmp = k_*spec3D.id(0);
         mat.second = Eigen::kroneckerProduct(tmp, spec1D.qDiff(2,0));

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
         mat.first = Eigen::kroneckerProduct(tmp, Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 2));

         // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void timeBlock(const BoussinesqBetaCylGTransport& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat k)
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

      // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
      mat.first = Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0));

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void boundaryBlock(const BoussinesqBetaCylGTransport& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      /// <b>Boundary operators for the transport equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$
      int pX = 0;
      int pZ = 0;

      MHDFloat cX;
      MHDFloat cZ;

      // Boundary condition for the streamfunction
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Set boundary condition prefactors
         cX = 0.0;
         cZ = 0.0;

      // Boundary condition for the vertical velocity
      } else if(fieldId.first == PhysicalNames::VELOCITYZ)
      {
         // Set boundary condition prefactors
         cX = 0.0;
         cZ = 0.0;

      // Boundary condition for the temperature
      } else if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         // Set boundary condition prefactors
         cX = 1.0;
         cZ = 0.0;

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for boundary operator!");
      }

      // Compute boundary block operator
      EquationToolsType::boundaryBlock(eq, FieldComponents::Spectral::SCALAR, mat, fieldId, pX, pZ, cX, cZ);
   }

}
}
