/** 
 * @file BoussinesqPerBetaCylGSystem.cpp
 * @brief Source of the implementation of the system of equations for the 3DQG beta model with periodic radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGSystem.hpp"

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

   void BoussinesqPerBetaCylGSystem::setCouplingInfo(CouplingInformation& rInfo, const SpectralFieldId eqId, const int nZ, const int modes)
   {
      /// - Streamfunction equation
      if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // General setup: first complex solver, complex solver, start from m = 0
         rInfo.setGeneral(0, true, 1);

         // 
         //  WARNING: the order is important as it determines the field index!
         //

         // Equation is coupled to streamfunction equation (self)
         rInfo.addImplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR, true);
         // Equation is coupled to vertical velocity equation
         rInfo.addImplicitField(PhysicalNames::VELOCITYZ,FieldComponents::Spectral::SCALAR, false);

         // Equation has explicit temperature
         //rInfo.addExplicitField(PhysicalNames::TEMPERATURE,FieldComponents::Spectral::SCALAR);

         // Set sizes of blocks and matrices
         ArrayI blockNs(modes);
         blockNs.setConstant(nZ);
         ArrayI rhsCols(modes);
         rhsCols.setConstant(1);
         rInfo.setSizes(modes, blockNs, rhsCols); 

      /// - Vertical velocity equation
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // General setup: first complex solver, complex solver, start from m = 0
         rInfo.setGeneral(0, true, 1);

         // 
         //  WARNING: the order is important as it determines the field index!
         //

         // Equation is coupled to streamfunction equation
         rInfo.addImplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR, false);
         // Equation is coupled to vertical velocity equation (self)
         rInfo.addImplicitField(PhysicalNames::VELOCITYZ,FieldComponents::Spectral::SCALAR, true);

         // Set sizes of blocks and matrices
         ArrayI blockNs(modes);
         blockNs.setConstant(nZ);
         ArrayI rhsCols(modes);
         rhsCols.setConstant(1);
         rInfo.setSizes(modes, blockNs, rhsCols); 

      /// - Transport equation
      } else if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         // General setup: first real solver, real solver, start from m = 0
         rInfo.setGeneral(0, false, 1);

         // 
         //  WARNING: the order is important
         //

         // Equation is coupled to temperature equation
         rInfo.addImplicitField(PhysicalNames::TEMPERATURE,FieldComponents::Spectral::SCALAR, true);

         // Equation has explicit temperature
         //rInfo.addExplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR);

         // Set sizes of blocks and matrices
         ArrayI blockNs(modes);
         blockNs.setConstant(nZ);
         ArrayI rhsCols(modes);
         rhsCols.setConstant(1);
         rInfo.setSizes(modes, blockNs, rhsCols); 

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID to set coupling information!");
      }
   }

   void BoussinesqPerBetaCylGSystem::quasiInverse(SparseMatrix& mat, const SpectralFieldId eqId, const int nZ)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nZ);

      /// - Streamfunction equation: \f$ D_Z^{-1} \f$
      if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         mat = spec1D.qDiff(1,0);

      /// - Vertical velocity equation: \f$ D_Z^{-1} \f$
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         mat = spec1D.qDiff(1,0);

      /// - Transport equation: \f$ I_Z \f$
      } else if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
         mat = spec1D.id(0);

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID for quasi-inverse operator!");
      }

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void BoussinesqPerBetaCylGSystem::timeBlock(DecoupledZSparse& mat, FieldComponents::Spectral::Id compId, const SpectralFieldId eqId, const int nZ, const MHDFloat kX, const MHDFloat kY, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nZ);

      // Initialise output matrices
      mat.first.resize(nZ,nZ);
      mat.second.resize(nZ,nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat kX_ = kX/2.;
      MHDFloat kY_ = kY/2.;

      /// - Streamfunction equation: \f$D_Z^{-1}\f$
      if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
         mat.first = Spectral::PeriodicOperator::laplacian2D(kX_, kY_)*spec1D.qDiff(1,0);

      /// - Vertical velocity equation: \f$D_Z^{-1}\f$
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         // Set time matrices (kronecker(A,B,out) => out = A(i,j)*B)
         mat.first = spec1D.qDiff(1,0);

      /// - Transport equation: \f$I_Z\f$
      } else if(eqId.first == PhysicalNames::TEMPERATURE)
      {
         // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
         mat.first = spec1D.id(0);

      // Unknown equation
      } else
      {
         throw Exception("Unknown equation ID for time operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void BoussinesqPerBetaCylGSystem::linearBlock(DecoupledZSparse& mat, FieldComponents::Spectral::Id compId, const SpectralFieldId eqId, const SpectralFieldId fieldId, const int nZ, const MHDFloat kX, const MHDFloat kY, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nZ);

      // Initialise output matrices
      mat.first.resize(nZ,nZ);
      mat.second.resize(nZ,nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat kX_ = kX/2.;
      MHDFloat kY_ = kY/2.;

      /// <b>Linear operators for the streamfunction equation</b>
      if(eqId.first == PhysicalNames::STREAMFUNCTION)
      {
         /// - Streamfunction : \f$ D_Z^{-1} \f$
         if(fieldId.first == PhysicalNames::STREAMFUNCTION)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.first = Spectral::PeriodicOperator::bilaplacian2D(kX_, kY_)*spec1D.qDiff(1,0);

         /// - Vertical velocity : \f$ I_z^{-1} \f$
         } else if(fieldId.first == PhysicalNames::VELOCITYZ)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.first = spec1D.id(1);

         /// - Temperature : \f$ i \frac{k}{2}\frac{1}{16}\frac{Ra}{Pr}D_Z^{-1} \f$
         } else if(fieldId.first == PhysicalNames::TEMPERATURE)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.second = kY_*(1.0/16.)*(Ra/Pr)*spec1D.qDiff(1,0);

         // Unknown field
         } else
         {
            throw Exception("Unknown field ID for linear operator!");
         }

      /// <b>Linear operators for the vertical velocity equation</b>
      } else if(eqId.first == PhysicalNames::VELOCITYZ)
      {
         /// - Streamfunction : \f$ -\frac{1}{\Gamma^2} I_Z^{-1} \f$
         if(fieldId.first == PhysicalNames::STREAMFUNCTION)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.first = -(1./(Gamma*Gamma))*spec1D.id(1);

         /// - Vertical velocity : \f$ D_Z^{-1} \f$
         } else if(fieldId.first == PhysicalNames::VELOCITYZ)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.first = Spectral::PeriodicOperator::laplacian2D(kX_, kY_)*spec1D.qDiff(1,0);

         /// - Temperature : \f$ 0_Z \f$
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
         /// - Streamfunction : \f$ i \frac{k}{2} I_Z \f$
         if(fieldId.first == PhysicalNames::STREAMFUNCTION)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.second = kY_*spec1D.id(0);

         /// - Vertical velocity : \f$ 0_Z \f$
         } else if(fieldId.first == PhysicalNames::VELOCITYZ)
         {
            //
            // Nothing to do for an empty sparse matrix
            //

         /// - Temperature : \f$ \frac{1}{Pr} I_Z \f$
         } else if(fieldId.first == PhysicalNames::TEMPERATURE)
         {
            // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
            mat.first = (1./Pr)*Spectral::PeriodicOperator::laplacian2D(kX_, kY_)*spec1D.id(0);

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

   void BoussinesqPerBetaCylGSystem::boundaryBlock(DecoupledZSparse& mat, FieldComponents::Spectral::Id compId, const SpectralFieldId eqId, const SpectralFieldId fieldId, const SharedSimulationBoundary spBcIds, const int nZ, const MHDFloat kX, const MHDFloat kY, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi)
   {
      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nZ);

      // Create spectral boundary operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::BcType bound1D(nZ);

      // Initialise output matrices
      mat.first.resize(nZ,nZ);
      mat.second.resize(nZ,nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat kY_ = kY/2.;

      // Create storage for the X boundary operators
      DecoupledZSparse tau;

      // Set boundary conditions on fieldId
      if(spBcIds->hasEquation(eqId) && spBcIds->hasField(eqId,fieldId))
      {
         // Set Z boundary conditions
         if(spBcIds->bcs(eqId,fieldId).count(Dimensions::Simulation::SIM1D) > 0)
         {
            tau = Spectral::BoundaryConditions::tauMatrix(bound1D, spBcIds->bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second);
            if(tau.first.nonZeros() > 0)
            {
               mat.first += tau.first;
            }
            if(tau.second.nonZeros() > 0)
            {
               tau.second *= kY_*std::tan((MathConstants::PI/180.)*chi)/Gamma;
               mat.second += tau.second;
            }
         }
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);

   }

   BoussinesqPerBetaCylGSystem::BoussinesqPerBetaCylGSystem()
   {
   }

   BoussinesqPerBetaCylGSystem::~BoussinesqPerBetaCylGSystem()
   {
   }
}
}
