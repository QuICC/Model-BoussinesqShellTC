/** \file TestTFTSystem.cpp
 *  \brief Source of the implementation of the system of equations for the TFT scheme tests
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
#include "Equations/Tests/TestTFTSystem.hpp"

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

   void TestTFTSystem::setCouplingInfo(CouplingInformation& rInfo, const SpectralFieldId eqId, const int nX, const int nZ, const int nY)
   {
      /// - First scalar equation
      if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         // General setup: first real solver, real solver, start from m = 0
         rInfo.setGeneral(0, false, 0);

         // 
         //  WARNING: the order is important
         //

         // Add self coupling
         rInfo.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR, true);

         // Set sizes of blocks and matrices
         ArrayI blockNs(nY);
         blockNs.setConstant(nX*nZ);
         ArrayI rhsCols(nY);
         rhsCols.setConstant(1);
         rInfo.setSizes(nY, blockNs, rhsCols); 

      /// - Second scalar equation
      } else if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // General setup: second real solver, real solver, start from m = 0
         rInfo.setGeneral(1, false, 0);

         // 
         //  WARNING: the order is important
         //

         // Add self coupling
         rInfo.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR, true);

         // Set sizes of blocks and matrices
         ArrayI blockNs(nY);
         blockNs.setConstant(nX*nZ);
         ArrayI rhsCols(nY);
         rhsCols.setConstant(1);
         rInfo.setSizes(nY, blockNs, rhsCols); 

      /// - Third scalar equation
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // General setup: third real solver, real solver, start from m = 0
         rInfo.setGeneral(2, false, 0);

         // 
         //  WARNING: the order is important
         //

         // Add self coupling
         rInfo.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR, true);

         // Set sizes of blocks and matrices
         ArrayI blockNs(nY);
         blockNs.setConstant(nX*nZ);
         ArrayI rhsCols(nY);
         rhsCols.setConstant(1);
         rInfo.setSizes(nY, blockNs, rhsCols); 

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID to set the coupling information!");
      }
   }

   void TestTFTSystem::quasiInverse(SparseMatrix& mat, const SpectralFieldId eqId, const int nX, const int nZ)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      /// - First scalar equation: \f$ \left( D_x^{-2} \otimes D_Z^{-2}\right) \f$
      if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(2,0), spec1D.qDiff(2,0), mat);

      /// - Second scalar equation: \f$ \left( D_x^{-2} \otimes D_Z^{-2}\right) \f$
      } else if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(2,0), spec1D.qDiff(2,0), mat);

      /// - Third scalar equation: \f$ \left( D_x^{-2} \otimes D_Z^{-2}\right) \f$
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(2,0), spec1D.qDiff(2,0), mat);

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID for quasi-inverse operator!");
      }

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void TestTFTSystem::timeBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const int nX, const int nZ, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Initialise output matrices
      mat.first.resize(nX*nZ,nX*nZ);
      mat.second.resize(nX*nZ,nX*nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      /// - First scalar equation: \f$\left(D_x^{-2} \otimes D_Z^{-2}\right)\f$
      if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(2,0), spec1D.qDiff(2,0), mat.first);

      /// - Second scalar equation: \f$\left(D_x^{-2} \otimes D_Z^{-2}\right)\f$
      } else if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(2,0), spec1D.qDiff(2,0), mat.first);

      /// - Third scalar equation: \f$\left(D_x^{-2} \otimes D_Z^{-2}\right)\f$
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(2,0), spec1D.qDiff(2,0), mat.first);

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID for time operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void TestTFTSystem::linearBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const int nX, const int nZ, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Initialise output matrices
      mat.first.resize(nX*nZ,nX*nZ);
      mat.second.resize(nX*nZ,nX*nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      /// <b>Linear operators for the first scalar equation</b>
      if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         /// - Streamfunction : \f$ \left(0_x \otimes 0_Z\right) \f$
         if(fieldId.first == PhysicalNames::STREAMFUNCTION)
         {
            //
            // Nothing to do for an empty sparse matrix
            //

         /// - Vertical velocity : \f$ \left(0_x \otimes 0_Z\right) \f$
         } else if(fieldId.first == PhysicalNames::VELOCITYZ)
         {
            //
            // Nothing to do for an empty sparse matrix
            //

         /// - Temperature : \f$ \left(D_x^{-2} \otimes D_Z^{-2}\right)\nabla^{2} \f$
         } else if(fieldId.first == PhysicalNames::TEMPERATURE)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.first = Spectral::PeriodicOperator::qLaplacian3D(spec1D, spec3D, k_, 2, 2);

         // Unknown field
         } else
         {
            throw Exception("Unknown field ID for linear operator!");
         }

      /// <b>Linear operators for the second scalar equation</b>
      } else if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         /// - Streamfunction : \f$ \left(0_x \otimes 0_Z\right) \f$
         if(fieldId.first == PhysicalNames::STREAMFUNCTION)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.first = Spectral::PeriodicOperator::qLaplacian3D(spec1D, spec3D, k_, 2, 2);

         /// - Vertical velocity : \f$ \left(0_x \otimes 0_Z\right) \f$
         } else if(fieldId.first == PhysicalNames::VELOCITYZ)
         {
            //
            // Nothing to do for an empty sparse matrix
            //

         /// - Temperature : \f$ \left(D_x^{-2} \otimes D_Z^{-2}\right)\nabla^{2} \f$
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

      /// <b>Linear operators for the first scalar equation</b>
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         /// - Streamfunction : \f$ \left(0_x \otimes 0_Z\right) \f$
         if(fieldId.first == PhysicalNames::STREAMFUNCTION)
         {
            //
            // Nothing to do for an empty sparse matrix
            //

         /// - Vertical velocity : \f$ \left(0_x \otimes 0_Z\right) \f$
         } else if(fieldId.first == PhysicalNames::VELOCITYZ)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.first = Spectral::PeriodicOperator::qLaplacian3D(spec1D, spec3D, k_, 2, 2);

         /// - Temperature : \f$ \left(D_x^{-2} \otimes D_Z^{-2}\right)\nabla^{2} \f$
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

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID for linear operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);

   }

   void TestTFTSystem::boundaryBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const SharedSimulationBoundary spBcIds, const int nX, const int nZ, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
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

      /// <b>Boundary operators for the first scalar equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$
      if(eqId.first == PhysicalNames::TEMPERATURE && spBcIds->hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.id(0);
         // Set Z boundary quasi-inverse
         q3D = spec3D.shiftId(2);

      /// <b>Boundary operators for the transport equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$
      } else if(eqId.first == PhysicalNames::STREAMFUNCTION && spBcIds->hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.id(0);
         // Set Z boundary quasi-inverse
         q3D = spec3D.shiftId(2);

      /// <b>Boundary operators for the transport equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$
      } else if(eqId.first == PhysicalNames::VELOCITYZ && spBcIds->hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.id(0);
         // Set Z boundary quasi-inverse
         q3D = spec3D.shiftId(2);

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
         // Set first dimension boundary conditions
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

         // Set third dimension boundary conditions
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

   TestTFTSystem::TestTFTSystem()
   {
   }

   TestTFTSystem::~TestTFTSystem()
   {
   }
}
}
