/** \file Beta3DQGStreamfunction.cpp
 *  \brief Source of the implementation of the streamfunction equation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGStreamfunction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   Beta3DQGStreamfunction::Beta3DQGStreamfunction(SharedIEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      /// \mhdBug Boundary conditions are not implemented correctly
      // Equation is always complex due to the sloping boundary condition
      this->setComplex(true);

      // Set the variable requirements
      this->setRequirements();
//
//         // Add boundary conditions for 1D
//            // ... No penetration 
//         this->addBC(FieldComponents::Spectral::SCALAR, Dimensions::ONED, std::make_pair(BoundaryConditions::VALUE,BoundaryConditions::BOTH));
//            // ... No-Slip
//         this->addBC(FieldComponents::Spectral::SCALAR, Dimensions::ONED, std::make_pair(BoundaryConditions::FIRST_DERIVATIVE,BoundaryConditions::BOTH));
//            // ... Stress-free
//         //this->addBC(FieldComponents::Spectral::SCALAR, Dimensions::ONED, std::make_pair(BoundaryConditions::SECOND_DERIVATIVE,BoundaryConditions::BOTH));
//         // Add boundary conditions for 3D
//         this->addBC(FieldComponents::Spectral::SCALAR, Dimensions::THREED, std::make_pair(BoundaryConditions::BETA_SLOPE,BoundaryConditions::TOP));
//         // Add coupled boundary conditions for 3D
//         this->addCBC(FieldComponents::Spectral::SCALAR, Dimensions::THREED, std::make_pair(BoundaryConditions::VALUE,BoundaryConditions::TOP));
   }

   Beta3DQGStreamfunction::~Beta3DQGStreamfunction()
   {
   }

   void Beta3DQGStreamfunction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      Physical::StreamAdvection::set(rNLComp, this->unknown().dom(0).grad(), this->scalar(PhysicalNames::VORTICITYZ).dom(0).grad(), 1.0);
   }

   void Beta3DQGStreamfunction::computeLinear(Datatypes::SpectralScalarType& rRHS) const
   {
      // Loop over m
      MHDFloat m_;

      // Compute Ra/(16 Pr As)
      MHDFloat c = -this->eqParams().nd(NonDimensional::RAYLEIGH)/(16.*this->eqParams().nd(NonDimensional::PRANDTL));

      for(int m = 0; m < this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         m_ = static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m));

         rRHS.addSlice((MathConstants::cI*m_*c)*this->scalar(PhysicalNames::TEMPERATURE).dom(0).perturbation().slice(m), m);
      }
   }

   void Beta3DQGStreamfunction::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::STREAMFUNCTION);

      // Set streamfunction requirements
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, true, true));

      // Add vertical velocity requirements
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, false, true));

      // Add temperature requirements
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, false, false, true));
   }

   void Beta3DQGStreamfunction::setCoupling()
   {
      // Set field coupling to vertical velocity
      this->mCouplingInfo.addField(FieldComponents::Spectral::SCALAR, std::make_pair(PhysicalNames::VELOCITYZ,FieldComponents::Spectral::SCALAR));

      // Set internal coupling (R and Z are coupled)
      int nMat = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      ArrayI dim(nMat);
      dim.setConstant(1);
      this->mCouplingInfo.addInternal(FieldComponents::Spectral::SCALAR, nMat, dim);
   }

   void Beta3DQGStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarEquation::timestepOutput(id, storage, matIdx, start);

      // Get right wave number
      int m_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(matIdx);

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Transform::TRA1D>::Type spec1D(this->unknown().dom(0).perturbation().slice(matIdx).rows());

      // Compute the vertical vorticity: \f$\eta = \nabla^2\psi\f$
      this->rScalar(PhysicalNames::VORTICITYZ).rDom(0).rPerturbation().setSlice(Spectral::PeriodicOperator::laplacian2D(spec1D, m_, 0)*this->unknown().dom(0).perturbation().slice(matIdx), matIdx);
   }

   void Beta3DQGStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarEquation::timestepOutput(id, storage, matIdx, start);

      // Get right wave number
      int m_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(matIdx);

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Transform::TRA1D>::Type spec1D(this->unknown().dom(0).perturbation().slice(matIdx).rows());

      // Compute the vertical vorticity: \f$\eta = \nabla^2\psi\f$
      this->rScalar(PhysicalNames::VORTICITYZ).rDom(0).rPerturbation().setSlice(Spectral::PeriodicOperator::laplacian2D(spec1D, m_, 0)*this->unknown().dom(0).perturbation().slice(matIdx), matIdx);
   }
 
   void Beta3DQGStreamfunction::setSpectralMatrices(Spectral::SpectralSelector<Dimensions::Transform::TRA1D>::Type& spec1D, Spectral::SpectralSelector<Dimensions::Transform::TRA2D>::Type& spec2D, Spectral::SpectralSelector<Dimensions::Transform::TRA3D>::Type& spec3D)
   {
//      // Simplify notation
//      typedef std::map<FieldComponents::Spectral::Id, std::vector<DecoupledZSparse> >::iterator MapItType;
//      typedef std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::iterator MapNLItType;
//      std::pair<MapItType,bool> pos;
//      MapItType it;
//      MapNLItType itNL;
//
//      // Temporary storage
//      SparseMatrix   tmpA;
//      SparseMatrix   tmpB;
//      DecoupledZSparse   tau1D;
//      DecoupledZSparse   tau3D;
//
//      bool atTop = true;
//
//      int bcSign;
//      if(atTop)
//      {
//         bcSign = 1;
//      } else
//      {
//         bcSign = -1;
//      }
//
//      // Get number of boundary conditions
//      int nBC1D = bcSign*this->nBC(FieldComponents::Spectral::SCALAR,Dimensions::ONED);
//      int nBC2D = bcSign*this->nBC(FieldComponents::Spectral::SCALAR,Dimensions::TWOD);
//      int nBC3D = bcSign*this->nBC(FieldComponents::Spectral::SCALAR,Dimensions::THREED);
//
//      //////////////////////////////////////////////
//      // Initialise boundary condition matrices
//      pos = this->mBCMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set boundary condition matrices
//      tau1D = spec1D.tau(this->getBCs(FieldComponents::Spectral::SCALAR,Dimensions::ONED), atTop);
//      Eigen::kroneckerProduct(spec2D.id(),spec3D.id(nBC3D),tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, tau1D.first, it->second.back().first);
//      tau3D.second = std::tan((MathConstants::PI/180.)*this->eqParams().nd(NonDimensional::CHI))*spec3D.tau(this->getBCs(FieldComponents::Spectral::SCALAR,Dimensions::THREED), atTop).second;
//      Eigen::kroneckerProduct(spec2D.id(),tau3D.second,tmpA);
//      Eigen::kroneckerProduct(tmpA, spec1D.id(), it->second.back().second);
//
//      // Set nonlinear multiplication matrix
//      itNL = this->mNLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<SparseMatrix>())).first;
//      Eigen::kroneckerProduct(spec2D.id(nBC2D),spec3D.qDiff(1,0),tmpA);
//      itNL->second.push_back(SparseMatrix());
//      Eigen::kroneckerProduct(tmpA, spec1D.qDiff(4,0), itNL->second.back());
//
//      //////////////////////////////////////////////
//      // Initialise coupled boundary condition matrices
//      pos = this->mCBCMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set coupled boundary condition matrices
//      tau3D = spec3D.tau(this->getCBCs(FieldComponents::Spectral::SCALAR,Dimensions::THREED),atTop);
//      Eigen::kroneckerProduct(spec2D.id(),tau3D.first,tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, spec1D.id(), it->second.back().first);
//
//      //////////////////////////////////////////////
//      // Initialise time matrices
//      pos = this->mTMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set time matrices
//      Eigen::kroneckerProduct(spec2D.id(),spec3D.qDiff(1,0),tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, spec1D.qPerpLaplacian(4,1), it->second.back().first);
//
//      //////////////////////////////////////////////
//      // Initialise linear matrices
//      pos = this->mLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set linear matrices
//      Eigen::kroneckerProduct(spec2D.id(),spec3D.qDiff(1,0),tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, spec1D.qPerpLaplacian(4,2), it->second.back().first);
//
//      //////////////////////////////////////////////
//      // Initialise coupling matrices
//      pos = this->mCMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set coupling matrices
//      tmpB = (-1.0/this->eqParams().nd(NonDimensional::GAMMA))*spec3D.id(nBC3D);
//      Eigen::kroneckerProduct(spec2D.id(),tmpB,tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, spec1D.qDiff(4,0), it->second.back().first);
   }
}
}
