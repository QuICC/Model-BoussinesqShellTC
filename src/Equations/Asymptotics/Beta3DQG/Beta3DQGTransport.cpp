/** \file Beta3DQGTransport.cpp
 *  \brief Source of the implementation of the transport equation in the 3DQG beta model
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
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamHeatAdvection.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   Beta3DQGTransport::Beta3DQGTransport(SharedIEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      /// \mhdBug Boundary conditions are not implemented correctly
      // Set the variable requirements
      this->setRequirements();
//
//         // Add boundary conditions for 1D
//         this->addBC(FieldComponents::Spectral::SCALAR, Dimensions::ONED, std::make_pair(BoundaryConditions::VALUE,BoundaryConditions::BOTH));
   }

   Beta3DQGTransport::~Beta3DQGTransport()
   {
   }

   void Beta3DQGTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      Physical::StreamHeatAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);
   }

   void Beta3DQGTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature to requirements
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, true));

      // Add streamfunction to requirements
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

   void Beta3DQGTransport::setCoupling()
   {
      // Set field coupling
      // NONE 

      // Set internal coupling (X and Z are coupled)
      int nMat = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      ArrayI dim(nMat);
      dim.setConstant(1);
      this->mCouplingInfo.addInternal(FieldComponents::Spectral::SCALAR, nMat, dim);
   }

   void Beta3DQGTransport::setSpectralMatrices(Spectral::SpectralSelector<Dimensions::Transform::TRA1D>::Type& spec1D, Spectral::SpectralSelector<Dimensions::Transform::TRA2D>::Type& spec2D, Spectral::SpectralSelector<Dimensions::Transform::TRA3D>::Type& spec3D)
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
//      Eigen::kroneckerProduct(spec2D.id(),spec3D.id(),tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, tau1D.first, it->second.back().first);
//
//      // Set nonlinear multiplication matrix
//      itNL = this->mNLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<SparseMatrix>())).first;
//      Eigen::kroneckerProduct(spec2D.id(nBC2D),spec3D.id(nBC3D),tmpA);
//      itNL->second.push_back(SparseMatrix());
//      Eigen::kroneckerProduct(tmpA, spec1D.qDiff(2,0), itNL->second.back());
//
//      //////////////////////////////////////////////
//      // Initialise time matrices
//      pos = this->mTMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set time matrices
//      Eigen::kroneckerProduct(spec2D.id(),spec3D.id(),tmpA);
//      it->second.push_back(DecoupledZSparse());
//      Eigen::kroneckerProduct(tmpA, spec1D.qDiff(2,0), it->second.back().first);
//
//      //////////////////////////////////////////////
//      // Initialise linear matrices
//      pos = this->mLMatrices.insert(std::make_pair(FieldComponents::Spectral::SCALAR, std::vector<DecoupledZSparse>()));
//      it = pos.first;
//
//      // Set linear matrices
//      Eigen::kroneckerProduct(spec2D.id(),spec3D.id(),tmpA);
//      it->second.push_back(DecoupledZSparse());
//      tmpB = (1.0/this->eqParams().nd(NonDimensional::PRANDTL))*spec1D.qPerpLaplacian(2,1);
//      Eigen::kroneckerProduct(tmpA, tmpB, it->second.back().first);
//
//      //////////////////////////////////////////////
//      // Initialise coupling matrices
//      // Set coupling matrices
//      // NONE
   }
}
}
