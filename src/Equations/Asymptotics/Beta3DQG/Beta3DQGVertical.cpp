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

namespace GeoMHDiSCC {

namespace Equations {

   Beta3DQGVertical::Beta3DQGVertical(SharedIEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      /// \mhdBug Boundary conditions are not implemented correctly
      // Equation is always complex due to the sloping boundary condition
      this->setComplex(true);

      // Set the variable requirements
      this->setRequirements();
//
//         // Add boundary conditions for 1D
//            // .. No-Slip
//         this->addBC(FieldComponents::Spectral::SCALAR, Dimensions::ONED, std::make_pair(BoundaryConditions::VALUE,BoundaryConditions::BOTH));
//            // .. Stress-free
//         //this->addBC(FieldComponents::Spectral::SCALAR, Dimensions::ONED, std::make_pair(BoundaryConditions::FIRST_DERIVATIVE,BoundaryConditions::BOTH));
//         // Add boundary conditions for 3D
//         this->addBC(FieldComponents::Spectral::SCALAR, Dimensions::THREED, std::make_pair(BoundaryConditions::VALUE,BoundaryConditions::BOTTOM));
//         // Add coupled boundary conditions for 3D
//         this->addCBC(FieldComponents::Spectral::SCALAR, Dimensions::THREED, std::make_pair(BoundaryConditions::BETA_SLOPE,BoundaryConditions::BOTTOM));
   }

   Beta3DQGVertical::~Beta3DQGVertical()
   {
   }

   void Beta3DQGVertical::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);
   }

   void Beta3DQGVertical::setRequirements()
   {
      // Set vertical velocity as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Add vertical velocity requirements
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));

      // Add streamfunction requirements
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

   void Beta3DQGVertical::setCoupling()
   {
      // Set field coupling to vertical velocity
      this->mCouplingInfo.addField(FieldComponents::Spectral::SCALAR, std::make_pair(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR));

      // Set internal coupling (R and Z are coupled)
      int nMat = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      ArrayI dim(nMat);
      dim.setConstant(1);
      this->mCouplingInfo.addInternal(FieldComponents::Spectral::SCALAR, nMat, dim);
   }
 
   void Beta3DQGVertical::setSpectralMatrices(Spectral::SpectralSelector<Dimensions::Transform::TRA1D>::Type& spec1D, Spectral::SpectralSelector<Dimensions::Transform::TRA2D>::Type& spec2D, Spectral::SpectralSelector<Dimensions::Transform::TRA3D>::Type& spec3D)
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
//      tau3D = spec3D.tau(this->getBCs(FieldComponents::Spectral::SCALAR,Dimensions::THREED), atTop);
//      Eigen::kroneckerProduct(spec2D.id(),tau3D.first,tmpA);
//      Eigen::kroneckerProduct(tmpA, spec1D.id(), tmpB);
//      it->second.back().first += tmpB;
//
//      // Set nonlinear multiplication matrix
//      itNL = this->mNLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<SparseMatrix>())).first;
//      Eigen::kroneckerProduct(spec2D.id(nBC2D),spec3D.qDiff(1,0),tmpA);
//      itNL->second.push_back(SparseMatrix());
//      Eigen::kroneckerProduct(tmpA, spec1D.qDiff(2,0), itNL->second.back());
//
//      //////////////////////////////////////////////
//      // Initialise coupled boundary condition matrices
//      pos = this->mCBCMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set coupled boundary condition matrices
//      tau3D.second = std::tan((MathConstants::PI/180.)*this->eqParams().nd(NonDimensional::CHI))*spec3D.tau(this->getCBCs(FieldComponents::Spectral::SCALAR,Dimensions::THREED), atTop).second;
//      Eigen::kroneckerProduct(spec2D.id(),tau3D.second,tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, spec1D.id(), it->second.back().second);
//
//      //////////////////////////////////////////////
//      // Initialise time matrices
//      pos = this->mTMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set time matrices
//      Eigen::kroneckerProduct(spec2D.id(),spec3D.qDiff(1,0),tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, spec1D.qDiff(2,0), it->second.back().first);
//
//      //////////////////////////////////////////////
//      // Initialise linear matrices
//      pos = this->mLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set linear matrices
//      Eigen::kroneckerProduct(spec2D.id(),spec3D.qDiff(1,0),tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, spec1D.qPerpLaplacian(2,1), it->second.back().first);
//
//      //////////////////////////////////////////////
//      // Initialise coupling matrices
//      pos = this->mCMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set coupling matrices
//      tmpB = (1.0/this->eqParams().nd(NonDimensional::GAMMA))*spec3D.id(nBC3D);
//      Eigen::kroneckerProduct(spec2D.id(),tmpB,tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, spec1D.qDiff(2,0), it->second.back().first);
   }
}
}
