/** \file Beta3DQGVertical.cpp
 *  \brief Source of the implementation of the vertical velocity equation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGVertical.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   Beta3DQGVertical::Beta3DQGVertical(SharedIEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Beta3DQGVertical::~Beta3DQGVertical()
   {
   }

   void Beta3DQGVertical::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///
      //Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);
      //
      rNLComp.rData().setConstant(0.0);
   }

   void Beta3DQGVertical::setRequirements()
   {
      // Equation is always complex due to the sloping boundary condition
      this->setComplex(true);

      // Set vertical velocity as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Add vertical velocity requirements
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));

      // Add streamfunction requirements
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

   void Beta3DQGVertical::setCoupling()
   {
      // Set the timestep starting index (exclude m = 0) mode
      if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0) == 0)
      {
         this->setStartIndex(1);
      }

      // Set field coupling to vertical velocity
      this->mCouplingInfo.addField(FieldComponents::Spectral::SCALAR, std::make_pair(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR));

      // Set internal coupling (R and Z are coupled)
      int nMat = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      ArrayI dim(nMat);
      dim.setConstant(1);
      this->mCouplingInfo.addInternal(FieldComponents::Spectral::SCALAR, nMat, dim);
   }
 
   void Beta3DQGVertical::setSpectralMatrices(const IEvolutionEquation::BcEqMapType& bcIds, const std::map<PhysicalNames::Id, IEvolutionEquation::BcEqMapType>& cbcIds)
   {
      // Get local copy of a shared resolution
      SharedResolution  spRes = this->unknown().dom(0).spRes();

      // Simplify notation
      typedef std::map<FieldComponents::Spectral::Id, std::vector<DecoupledZSparse> >::iterator spectral_iterator;
      typedef std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::iterator nonlin_iterator;

      // Create map position object
      std::pair<spectral_iterator,bool> pos;
      // Create spectral map iterator
      spectral_iterator it;
      // Create nonlinear map iterator
      nonlin_iterator itNL;

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(1);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(1);

      // Create spectral boundary operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::BcType bound1D(1);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::BcType bound3D(1);

      // Create boundary conditions ID
      IEvolutionEquation::BcKeyType bc1D = std::make_pair(FieldComponents::Spectral::SCALAR, Dimensions::Simulation::SIM1D);
      IEvolutionEquation::BcKeyType bc3D = std::make_pair(FieldComponents::Spectral::SCALAR, Dimensions::Simulation::SIM3D);

      // Temporary storage
      SparseMatrix   tmpA;
      DecoupledZSparse   tau1D;
      DecoupledZSparse   tau3D;

      // Get the box scale
      MHDFloat boxScale = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);

      // Get third dimension size
      int dim3D = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();

      /// 
      /// The timestepper will buid the full implicit timestepping matrices as L - T, where L is the linear operator and T is the time derivative matrix.
      /// The  left-hand side following equation is setup here: 
      ///
      ///    \f$ \nabla_{\perp}^2 w - \partial_t w + \frac{1}{\Gamma}\partial_Z \psi = \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///

      // Loop over third dimension
      for(int k = 0; k < dim3D; k++)
      {
         // Get global index in third data dimension (second physical dimension)
         MHDFloat k_ = boxScale*0.5*static_cast<MHDFloat>(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)); 

         // Reset spectral operator 1D
         spec1D.reset(spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
         // Reset spectral operator 3D
         spec3D.reset(spRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL));
         // Reset spectral boundary operator 1D
         bound1D.reset(spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
         // Reset spectral boundary operator 3D
         bound3D.reset(spRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL));

         //////////////////////////////////////////////
         // Initialise boundary condition matrices
         pos = this->mBCMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
         it = pos.first;

         // Set boundary condition matrices (kronecker(A,B,out) => out = A(i,j)*A)
         it->second.push_back(DecoupledZSparse());
         tau1D = Spectral::BoundaryConditions::tauMatrix(bound1D, bcIds.find(bc1D)->second);
//         Eigen::kroneckerProduct(spec3D.qDiff(1,0), tau1D.first, it->second.back().first);
         Eigen::kroneckerProduct(spec3D.id(0), tau1D.first, it->second.back().first);

         tau3D = Spectral::BoundaryConditions::tauMatrix(bound3D, bcIds.find(bc3D)->second);
//         Eigen::kroneckerProduct(tau3D.first, spec1D.id(0), tmpA);
         Eigen::kroneckerProduct(tau3D.first, spec1D.qDiff(2,0), tmpA);
         it->second.back().first += tmpA;

         // Prune matrices for safety
         it->second.back().first.prune(1e-32);
         it->second.back().second.prune(1e-32);

         //////////////////////////////////////////////
         // Initialise coupled boundary condition matrices
         pos = this->mCBCMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
         it = pos.first;

         // Set coupled boundary condition matrices (kronecker(A,B,out) => out = A(i,j)*A)
         it->second.push_back(DecoupledZSparse());
         tau3D = Spectral::BoundaryConditions::tauMatrix(bound3D, cbcIds.find(PhysicalNames::STREAMFUNCTION)->second.find(bc3D)->second);
         tau3D.second *= k_*std::tan((MathConstants::PI/180.)*this->eqParams().nd(NonDimensional::CHI));
//         Eigen::kroneckerProduct(tau3D.second, spec1D.id(0), it->second.back().second);
         Eigen::kroneckerProduct(tau3D.second, spec1D.qDiff(2,0), it->second.back().second);

         // Prune matrices for safety
         it->second.back().first.prune(1e-32);
         it->second.back().second.prune(1e-32);

         //////////////////////////////////////////////
         // Initialise nonlinear multiplication matrix
         itNL = this->mNLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<SparseMatrix>())).first;
         itNL->second.push_back(SparseMatrix());

         // Set nonlinear multiplication matrix (kronecker(A,B,out) => out = A(i,j)*A)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(2,0), itNL->second.back());

         // Prune matrices for safety
         itNL->second.back().prune(1e-32);
         itNL->second.back().prune(1e-32);

         //////////////////////////////////////////////
         // Initialise time matrices
         pos = this->mTMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
         it = pos.first;

         // Set time matrices (kronecker(A,B,out) => out = A(i,j)*A)
         it->second.push_back(DecoupledZSparse());
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(2,0), it->second.back().first);

         // Prune matrices for safety
         it->second.back().first.prune(1e-32);
         it->second.back().second.prune(1e-32);

         //////////////////////////////////////////////
         // Initialise linear matrices
         pos = this->mLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
         it = pos.first;

         // Set linear matrices (kronecker(A,B,out) => out = A(i,j)*A)
         it->second.push_back(DecoupledZSparse());
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 2), it->second.back().first);

         // Prune matrices for safety
         it->second.back().first.prune(1e-32);
         it->second.back().second.prune(1e-32);

         //////////////////////////////////////////////
         // Initialise coupling matrices
         pos = this->mCMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
         it = pos.first;

         // Set coupling matrices (kronecker(A,B,out) => out = A(i,j)*A)
         it->second.push_back(DecoupledZSparse());
         tmpA = (1.0/this->eqParams().nd(NonDimensional::GAMMA))*spec3D.id(1);
         Eigen::kroneckerProduct(tmpA, spec1D.qDiff(2,0), it->second.back().first);

         // Prune matrices for safety
         it->second.back().first.prune(1e-32);
         it->second.back().second.prune(1e-32);
      }
   }
}
}
