/** 
 * @file BoussinesqBetaCylGStreamfunction.cpp
 * @brief Source of the implementation of the streamfunction equation in the 3DQG beta model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "SpectralOperators/Tools/SpectralBoxTools.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

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

   void BoussinesqBetaCylGStreamfunction::setCoupling()
   {
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

      // Set mininal matrix coupling
      int nMat = 0;
      ArrayI blockNs;
      ArrayI rhsCols;
      EigenSelector::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

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

   DecoupledZSparse BoussinesqBetaCylGStreamfunction::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary) const
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

   void BoussinesqBetaCylGStreamfunction::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      quasiInverseBlock(*this, compId, mat);
   }

   void BoussinesqBetaCylGStreamfunction::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      linearBlock(*this, compId, mat, fieldId, eigs, false);
   }

   void quasiInverseBlock(const BoussinesqBetaCylGStreamfunction& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nX);
      Spectral::OperatorSelector<Dimensions::Simulation::SIM3D>::Type spec3D(nZ);

      EigenSelector::KRSum blocks;
      EigenSelector::KRProduct kProduct(SparseMatrix(nX*nZ,nX*nZ),SparseMatrix(nX*nZ,nX*nZ));

      /// - Streamfunction equation: \f$ \left(D_x^{-4} \otimes D_Z^{-1}\right) \f$
      std::tr1::get<0>(kProduct) = spec1D.qDiff(4,0);
      std::tr1::get<1>(kProduct) = spec3D.qDiff(1,0);
      blocks.push_back(kProduct);

      EigenSelector::computeKSum(mat, blocks);
   }

   void linearBlock(const BoussinesqBetaCylGStreamfunction& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
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

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = eq.eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = eq.eqParams().nd(NonDimensional::PRANDTL);

      /// - Streamfunction : \f$ \left(D_x^{-4}\nabla_\perp^{4} \otimes D_Z^{-1}\right) \f$
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         std::tr1::get<0>(kProduct).real() = Spectral::BoxTools::qBilaplacian2D(spec1D, k_, 4);
         std::tr1::get<1>(kProduct).real() = spec3D.qDiff(1,0);
         blocks.push_back(kProduct);

      /// - Vertical velocity : \f$ \left(D_x^{-4} \otimes I_z^{-1}\right) \f$
      } else if(fieldId.first == PhysicalNames::VELOCITYZ)
      {
         std::tr1::get<0>(kProduct).real() = spec1D.qDiff(4,0);
         std::tr1::get<1>(kProduct).real() = spec3D.id(1);
         blocks.push_back(kProduct);

         /// - Temperature : \f$ i \frac{k}{2}\frac{1}{16}\frac{Ra}{Pr}\left( D_x^{-4} \otimes D_Z^{-1}\right) \f$
      } else if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         std::tr1::get<0>(kProduct).imag() = spec1D.qDiff(4,0);
         std::tr1::get<1>(kProduct).imag() = k_*(1.0/16.)*(Ra/Pr)*spec3D.qDiff(1,0);
         blocks.push_back(kProduct);

         // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }

      EigenSelector::constrainBlock(eq, compId, mat, fieldId, blocks, eigs, hasBoundary);
   }

   void timeBlock(const BoussinesqBetaCylGStreamfunction& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
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

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         std::tr1::get<0>(kProduct).real() = Spectral::BoxTools::qLaplacian2D(spec1D, k_, 4);
         std::tr1::get<1>(kProduct).real() = spec3D.qDiff(1,0);
         blocks.push_back(kProduct);
      } else
      {
         throw Exception("Multiple field in time integration not implemented yet!");
      }

      EigenSelector::constrainBlock(eq, compId, mat, fieldId, blocks, eigs, hasBoundary);
   }

   void boundaryBlock(const BoussinesqBetaCylGStreamfunction& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx)
   {
      assert(eigs.size() == 1);
      MHDFloat k = eigs.at(0);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Get Physical parameters Gamma, chi
      MHDFloat Gamma = eq.eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = eq.eqParams().nd(NonDimensional::CHI);

      // Boundary condition prefactors for the streamfunction
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         coeffs.push_back(1.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

         coeffs.push_back(k_*std::tan((Math::PI/180.)*chi)/Gamma);
         bcIdx.push_back(Boundary::BCIndex(k));

      // Boundary condition prefactors for the vertical velocity
      } else if(fieldId.first == PhysicalNames::VELOCITYZ)
      {
         coeffs.push_back(1.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

         coeffs.push_back(1.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

      // Boundary condition prefactors for the temperature
      } else if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         coeffs.push_back(0.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

         coeffs.push_back(0.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for boundary operator!");
      }
   }

}
}
