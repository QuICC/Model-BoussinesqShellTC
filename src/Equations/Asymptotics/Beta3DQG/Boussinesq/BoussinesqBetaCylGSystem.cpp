/** \file BoussinesqBetaCylGSystem.cpp
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGSystem.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/Dimensions.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Simulation/SimulationBoundary.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   void BoussinesqBetaCylGSystem::quasiInverse(SparseMatrix& mat, const SpectralFieldId eqId, const int nX, const int nZ)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      /// - Streamfunction equation: \f$ \left(D_x^{-4} \otimes D_Z^{-1}\right) \f$
      if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(4,0), mat);

      /// - Vertical velocity equation: \f$ \left(D_x^{-2} \otimes D_Z^{-1}\right) \f$
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(2,0), mat);

      /// - Transport equation: \f$ \left( D_x^{-2} \otimes I_Z\right) \f$
      } else if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0), mat);

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID for quasi-inverse operator!");
      }

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void BoussinesqBetaCylGSystem::boundaryBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const SharedSimulationBoundary spBcIds, const int nX, const int nZ, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
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

      /// <b>Boundary operators for the streamfunction equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$ and \f$\left(D_x^{-4} \otimes BC_Z\right)\f$
      if(eqId.first == PhysicalNames::STREAMFUNCTION && spBcIds->hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.shiftId(4);

         // Set Z boundary quasi-inverse
         q3D = spec3D.id(0);

      /// <b>Boundary operators for the vertical velocity equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$ and \f$\left(D_x^{-2} \otimes BC_Z\right)\f$
      } else if(eqId.first == PhysicalNames::VELOCITYZ && spBcIds->hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.shiftId(2);
         // Set Z boundary quasi-inverse
         q3D = spec3D.id(0);

      /// <b>Boundary operators for the transport equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$
      } else if(eqId.first == PhysicalNames::TEMPERATURE && spBcIds->hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.id(0);
         // Set Z boundary quasi-inverse
         q3D = spec3D.id(0);

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID or missing condition for boundary operator!");
      }

      // Create storage for the X boundary operators
      DecoupledZSparse tau;

      // Set boundary conditions on fieldId
      if(spBcIds->hasField(eqId,fieldId))
      {
         // Set X boundary conditions
         if(spBcIds->bcs(eqId,fieldId).count(Dimensions::Simulation::SIM1D) > 0)
         {
            tau = Spectral::BoundaryConditions::tauMatrix(bound1D, spBcIds->bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second);
            if(tau.first.nonZeros() > 0)
            {
               Eigen::kroneckerProduct(q3D, tau.first, mat.first);
            }

            if(tau.second.nonZeros() > 0)
            {
               Eigen::kroneckerProduct(q3D, tau.second, mat.second);
            }
         }

         // Set Z boundary conditions
         if(spBcIds->bcs(eqId,fieldId).count(Dimensions::Simulation::SIM3D) > 0)
         {
            tau = Spectral::BoundaryConditions::tauMatrix(bound3D, spBcIds->bcs(eqId,fieldId).find(Dimensions::Simulation::SIM3D)->second);
            if(tau.first.nonZeros() > 0)
            {
               SparseMatrix tmp;
               Eigen::kroneckerProduct(tau.first, q1D, tmp);
               mat.first += tmp;
            }
            if(tau.second.nonZeros() > 0)
            {
               tau.second *= k_*std::tan((MathConstants::PI/180.)*chi)/Gamma;
               SparseMatrix tmp;
               Eigen::kroneckerProduct(tau.second, q1D, tmp);
               mat.second += tmp;
            }
         }
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }
}
}
