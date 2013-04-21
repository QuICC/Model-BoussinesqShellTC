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
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/Dimensions.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Simulation/SimulationBoundary.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   void Beta3DQGSystem::setCouplingInfo(CouplingInformation& rInfo, const SpectralFieldId eqId, const int nx, const int nz, const int ny)
   {
      /// - Streamfunction equation
      if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Generat setup: first complex solver, complex solver, start from m = 0
         rInfo.setGeneral(0, true, 0);

         // 
         //  WARNING: the order is important as it determines the field index!
         //

         // Equation is coupled to streamfunction equation (self)
         rInfo.addImplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR, true);
         // Equation is coupled to vertical velocity equation
         rInfo.addImplicitField(PhysicalNames::VELOCITYZ,FieldComponents::Spectral::SCALAR, false);

         // Equation has explicit temperature
         rInfo.addExplicitField(PhysicalNames::TEMPERATURE,FieldComponents::Spectral::SCALAR);

         // Set sizes of blocks and matrices
         ArrayI blockNs(ny);
         blockNs.setConstant(nx*nz);
         ArrayI rhsCols(ny);
         rhsCols.setConstant(1);
         rInfo.setSizes(ny, blockNs, rhsCols); 

      /// - Vertical velocity equation
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // Generat setup: first complex solver, complex solver, start from m = 0
         rInfo.setGeneral(0, true, 0);

         // 
         //  WARNING: the order is important as it determines the field index!
         //

         // Equation is coupled to streamfunction equation
         rInfo.addImplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR, false);
         // Equation is coupled to vertical velocity equation (self)
         rInfo.addImplicitField(PhysicalNames::VELOCITYZ,FieldComponents::Spectral::SCALAR, true);

         // Set sizes of blocks and matrices
         ArrayI blockNs(ny);
         blockNs.setConstant(nx*nz);
         ArrayI rhsCols(ny);
         rhsCols.setConstant(1);
         rInfo.setSizes(ny, blockNs, rhsCols); 

      /// - Transport equation
      } else if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         // Generat setup: first real solver, real solver, start from m = 0
         rInfo.setGeneral(0, false, 0);

         // 
         //  WARNING: the order is important
         //

         // Equation is coupled to temperature equation
         rInfo.addImplicitField(PhysicalNames::TEMPERATURE,FieldComponents::Spectral::SCALAR, true);

         // Equation has explicit temperature
         rInfo.addExplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR);

         // Set sizes of blocks and matrices
         ArrayI blockNs(ny);
         blockNs.setConstant(nx*nz);
         ArrayI rhsCols(ny);
         rhsCols.setConstant(1);
         rInfo.setSizes(ny, blockNs, rhsCols); 

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID for quasi-inverse operator!");
      }
   }

   void Beta3DQGSystem::quasiInverse(SparseMatrix& mat, const SpectralFieldId eqId, const int nx, const int nz)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nx);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nz);

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

   void Beta3DQGSystem::timeBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const int nx, const int nz, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nx);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nz);

      // Initialise output matrices
      mat.first.resize(nx*nz,nx*nz);
      mat.second.resize(nx*nz,nx*nz);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      /// - Streamfunction equation: \f$\left(D_x^{-4}\nabla_\perp^2 \otimes D_Z^{-1}\right)\f$
      if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 4), mat.first);

      /// - Vertical velocity equation: \f$\left(D_x^{-2} \otimes D_Z^{-1}\right)\f$
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // Set time matrices (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(2,0), mat.first);

      /// - Transport equation: \f$\left(D_x^{-2} \otimes I_Z\right)\f$
      } else if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.id(0), spec1D.qDiff(2,0), mat.first);

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID for time operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void Beta3DQGSystem::linearBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const int nx, const int nz, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nx);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nz);

      // Initialise output matrices
      mat.first.resize(nx*nz,nx*nz);
      mat.second.resize(nx*nz,nx*nz);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      /// <b>Linear operators for the streamfunction equation</b>
      if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
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

      /// <b>Linear operators for the vertical velocity equation</b>
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         /// - Streamfunction : \f$ -\frac{1}{\Gamma^2}\left(D_x^{-2} \otimes I_Z^{-1}\right) \f$
         if(fieldId.first == PhysicalNames::STREAMFUNCTION)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            SparseMatrix tmp = -(1./(Gamma*Gamma))*spec3D.id(1);
            Eigen::kroneckerProduct(tmp, spec1D.qDiff(2,0), mat.first);

         /// - Vertical velocity : \f$ \left(D_x^{-2}\nabla_\perp^{2}\otimes D_Z^{-1}\right) \f$
         } else if(fieldId.first == PhysicalNames::VELOCITYZ)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            Eigen::kroneckerProduct(spec3D.qDiff(1,0), Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 2), mat.first);

         /// - Temperature : \f$ \left(0_x \otimes 0_Z\right) \f$
         } else if(fieldId.first == PhysicalNames::TEMPERATURE)
         {
            //
            // Nothing to do for an empty sparse matrix
            //

         // Unknown field
         } else
         {
            throw Exception("Unknown field ID for linear operator!");
         }

      /// <b>Linear operators for the transport equation</b>
      } else if(eqId.first == PhysicalNames::TEMPERATURE)
      {
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

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID for linear operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);

   }

   void Beta3DQGSystem::boundaryBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const SharedSimulationBoundary spBcIds, const int nx, const int nz, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nx);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nz);

      // Create spectral boundary operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::BcType bound1D(nx);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::BcType bound3D(nz);

      // Initialise output matrices
      mat.first.resize(nx*nz,nx*nz);
      mat.second.resize(nx*nz,nx*nz);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Storage for the boundary quasi-inverses
      SparseMatrix q1D;
      SparseMatrix q3D;

      /// <b>Boundary operators for the streamfunction equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$ and \f$\left(D_x^{-4} \otimes BC_Z\right)\f$
      if(eqId.first == PhysicalNames::STREAMFUNCTION && spBcIds->hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.qDiff(4,0);
         // Set Z boundary quasi-inverse
         q3D = spec3D.id(0);

      /// <b>Boundary operators for the vertical velocity equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$ and \f$\left(D_x^{-2} \otimes BC_Z\right)\f$
      } else if(eqId.first == PhysicalNames::VELOCITYZ && spBcIds->hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.qDiff(2,0);
         // Set Z boundary quasi-inverse
         q3D = spec3D.id(0);

      /// <b>Boundary operators for the transport equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$
      } else if(eqId.first == PhysicalNames::TEMPERATURE && spBcIds->hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.qDiff(2,0);
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

   Beta3DQGSystem::Beta3DQGSystem()
   {
   }

   Beta3DQGSystem::~Beta3DQGSystem()
   {
   }
}
}
