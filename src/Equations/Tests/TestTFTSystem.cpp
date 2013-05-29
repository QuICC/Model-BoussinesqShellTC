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
         // Generat setup: first real solver, real solver, start from m = 0
         rInfo.setGeneral(0, false, 0);

         // 
         //  WARNING: the order is important
         //

         // Add self coupling
         rInfo.addImplicitField(PhysicalNames::TEMPERATURE,FieldComponents::Spectral::SCALAR, true);

         // Set sizes of blocks and matrices
         ArrayI blockNs(nY);
         blockNs.setConstant(nX*nZ);
         ArrayI rhsCols(nY);
         rhsCols.setConstant(1);
         rInfo.setSizes(nY, blockNs, rhsCols); 

      /// - Second scalar equation
      if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Generat setup: first real solver, real solver, start from m = 0
         rInfo.setGeneral(0, false, 0);

         // 
         //  WARNING: the order is important
         //

         // Add self coupling
         rInfo.addImplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR, true);

         // Set sizes of blocks and matrices
         ArrayI blockNs(nY);
         blockNs.setConstant(nX*nZ);
         ArrayI rhsCols(nY);
         rhsCols.setConstant(1);
         rInfo.setSizes(nY, blockNs, rhsCols); 

      /// - Third scalar equation
      if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // Generat setup: first real solver, real solver, start from m = 0
         rInfo.setGeneral(0, false, 0);

         // 
         //  WARNING: the order is important
         //

         // Add self coupling
         rInfo.addImplicitField(PhysicalNames::VELOCITYZ,FieldComponents::Spectral::SCALAR, true);

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

      /// - Transport equation: \f$ \left( D_x^{-2} \otimes D_Z^{-2}\right) \f$
      if(eqId.first == PhysicalNames::TEMPERATURE)
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

      /// - Transport equation: \f$\left(D_x^{-2} \otimes D_Z^{-2}\right)\f$
      if(eqId.first == PhysicalNames::TEMPERATURE)
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

      /// <b>Linear operators for the transport equation</b>
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
            Eigen::kroneckerProduct(spec3D.qDiff(2,0), spec1D.id(2), mat.first);
            SparseMatrix tmp;
            Eigen::kroneckerProduct(spec3D.qDiff(2,0), spec1D.qDiff(2,0), tmp);
            mat.first = mat.first - (k_*k_)*tmp;
            Eigen::kroneckerProduct(spec3D.id(2), spec1D.qDiff(2,0), tmp);
            mat.first = mat.first + tmp;

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

      /// <b>Boundary operators for the transport equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$
      if(eqId.first == PhysicalNames::TEMPERATURE && spBcIds->hasEquation(eqId))
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

   TestTFTSystem::TestTFTSystem()
   {
   }

   TestTFTSystem::~TestTFTSystem()
   {
   }
}
}
