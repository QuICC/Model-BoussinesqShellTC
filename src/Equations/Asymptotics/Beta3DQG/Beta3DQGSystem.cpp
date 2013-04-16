/** \file Beta3DQGSystem.cpp
 *  \brief Source of the implementation of the system of equations for the 3DQG beta model
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
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGSystem.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   Beta3DQGSystem::Beta3DQGSystem(SharedIEquationParameters spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Beta3DQGSystem::~Beta3DQGSystem()
   {
   }

   void Beta3DQGSystem::time() const
   {
      // Set time evolution operator of streamfunction equation
      if(eqId == PhysicalNames::STREAMFUNCTION)
      {
         // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 4), op.first);

         // Prune matrices for safety
         op.first.prune(1e-32);
         op.second.prune(1e-32);

      // Set time evolution operator of vertical velocity equation
      } else if(eqId == PhysicalNames::VELOCITYZ)
      {
         // Set time matrices (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(2,0), op.first);

         // Prune matrices for safety
         op.first.prune(1e-32);
         op.second.prune(1e-32);

      // Set time evolution operator of transport equation
      } else if(eqId == PhysicalNames::TEMPERATURE)
      {
         // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0), op.first);

         // Prune matrices for safety
         op.first.prune(1e-32);
         op.second.prune(1e-32);
      } else
      {
      }
   }

   void Beta3DQGSystem::quasiInverse(SparseMatrix& op, const PhysicalNames::Id eqId) const
   {
      // Set quasi-inverse operator of streamfunction equation
      if(eqId == PhysicalNames::STREAMFUNCTION)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(4,0), op);

         // Prune matrices for safety
         op.prune(1e-32);

      // Set quasi-inverse operator of vertical velocity equation
      } else if(eqId == PhysicalNames::VELOCITYZ)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(2,0), op);

         // Prune matrices for safety
         op.prune(1e-32);

      // Set quasi-inverse operator of transport equation
      } else if(eqId == PhysicalNames::TEMPERATURE)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0), op);

         // Prune matrices for safety
         op.prune(1e-32);
      } else
      {
      }
   }

   void Beta3DQGSystem::linearStreamfunction() const
   {
      if(fieldId == PhysicalNames::STREAMFUNCTION)
      {
         // Set streamfunction linear operator of the streamfunction equation (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), Spectral::PeriodicOperator::qBilaplacian2D(spec1D, k_, 4), op.first);

         // Prune matrices for safety
         op.first.prune(1e-32);
         op.second.prune(1e-32);

      } else if(fieldId == PhysicalNames::VELOCITYZ)
      {
         // Set vertical velocity linear operator of the streamfunction equation (kronecker(A,B,out) => out = A(i,j)*B)
         SparseMatrix tmp = (-1.0/eqParams.nd(NonDimensional::GAMMA))*spec3D.id(1);
         Eigen::kroneckerProduct(tmpA, spec1D.qDiff(4,0), op.first);

         // Prune matrices for safety
         op.first.prune(1e-32);
         op.second.prune(1e-32);

      } else if(fieldId == PhysicalNames::TEMPERATURE)
      {
         // Set temperature linear operator of the streamfunction equation (kronecker(A,B,out) => out = A(i,j)*B)
         SparseMatrix tmp = eqParams.nd(NonDimensional::RAYLEIGH)*spec3D.id(1);
         Eigen::kroneckerProduct(tmpA, spec1D.id(4), op.first);

         // Prune matrices for safety
         op.first.prune(1e-32);
         op.second.prune(1e-32);

      } else
      {
      }
   }

   void Beta3DQGSystem::linearVelocityZ() const
   {
      if(fieldId == PhysicalNames::STREAMFUNCTION)
      {
      } else if(fieldId == PhysicalNames::VELOCITYZ)
      {
      } else if(fieldId == PhysicalNames::TEMPERATURE)
      {
      } else
      {
      }
   }

   void Beta3DQGSystem::linearTemperature() const
   {
      if(fieldId == PhysicalNames::STREAMFUNCTION)
      {
      } else if(fieldId == PhysicalNames::VELOCITYZ)
      {
      } else if(fieldId == PhysicalNames::TEMPERATURE)
      {
      } else
      {
      }
   }

   void Beta3DQGSystem::boundaryStreamfunction() const
   {
      if(fieldId == PhysicalNames::STREAMFUNCTION)
      {
         DecoupledZSparse tau1D = Spectral::BoundaryConditions::tauMatrix(bound1D, bcIds.find(bc1D)->second);
         Eigen::kroneckerProduct(spec3D.id(0), tau1D.first, op.first);
         tau3D = Spectral::BoundaryConditions::tauMatrix(bound3D, bcIds.find(bc3D)->second);
         tau3D.second *= k_*std::tan((MathConstants::PI/180.)*this->eqParams().nd(NonDimensional::CHI));
         Eigen::kroneckerProduct(tau3D.second, spec1D.qDiff(4,0), op.second);

         // Prune matrices for safety
         op.first.prune(1e-32);
         op.second.prune(1e-32);
      } else if(fieldId == PhysicalNames::VELOCITYZ)
      {
         DecoupledZSparse tau3D = Spectral::BoundaryConditions::tauMatrix(bound3D, cbcIds.find(PhysicalNames::VELOCITYZ)->second.find(bc3D)->second);
         Eigen::kroneckerProduct(tau3D.first, spec1D.qDiff(4,0), it->second.back().first);

         // Prune matrices for safety
         op.first.prune(1e-32);
         op.second.prune(1e-32);

      } else if(fieldId == PhysicalNames::TEMPERATURE)
      {

         // Prune matrices for safety
         op.first.prune(1e-32);
         op.second.prune(1e-32);

      } else
      {
      }
   }

   void Beta3DQGSystem::boundaryVelocityZ() const
   {
      if(fieldId == PhysicalNames::STREAMFUNCTION)
      {
      } else if(fieldId == PhysicalNames::VELOCITYZ)
      {
      } else if(fieldId == PhysicalNames::TEMPERATURE)
      {
      } else
      {
      }
   }

   void Beta3DQGSystem::boundaryTemperature() const
   {
      if(fieldId == PhysicalNames::STREAMFUNCTION)
      {
      } else if(fieldId == PhysicalNames::VELOCITYZ)
      {
      } else if(fieldId == PhysicalNames::TEMPERATURE)
      {
      } else
      {
      }
   }
}
}
