/** \file Scalar1DEigenTools.cpp
 *  \brief Source of the tools for schemes with a single eigen direction
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
#include "Equations/Tools/Scalar1DEigenTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   void boundaryBlock1DPeriodic(const Scalar1DEigenTools& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int pX, const int pZ, const MHDFloat cX, const MHDFloat cZ)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create equation ID
      SpectralFieldId eqId = std::make_pair(eq.name(), FieldComponents::Spectral::SCALAR);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Create spectral boundary operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::BcType bound1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::BcType bound3D(nZ);

      // Initialise output matrices
      mat.first.resize(nX*nZ,nX*nZ);
      mat.second.resize(nX*nZ,nX*nZ);

      // Storage for the boundary quasi-inverses
      SparseMatrix q1D;
      SparseMatrix q3D;

      // Set boundary "operators"
      if(eq.bcIds().hasEquation(eqId))
      {
         // Set X boundary quasi-inverse
         q1D = spec1D.shiftId(pX);

         // Set Z boundary quasi-inverse
         q3D = spec3D.id(pZ);
      // Unknown equation
      } else
      {
         throw Exception("Missing condition(s) for boundary operator!");
      }

      // Temporary storage for the boundary operators
      DecoupledZSparse tau;

      // Set boundary conditions on fieldId
      if(eq.bcIds().hasField(eqId,fieldId))
      {
         // Set X boundary conditions
         if(eq.bcIds().bcs(eqId,fieldId).count(Dimensions::Simulation::SIM1D) > 0)
         {
            tau = Spectral::BoundaryConditions::tauMatrix(bound1D, eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second);
            if(tau.first.nonZeros() > 0)
            {
               if(cX != 1.0)
               {
                  tau.first *= cX;
               }
               Eigen::kroneckerProduct(q3D, tau.first, mat.first);
            }

            if(tau.second.nonZeros() > 0)
            {
               if(cX != 1.0)
               {
                  tau.second *= cX;
               }
               Eigen::kroneckerProduct(q3D, tau.second, mat.second);
            }
         }

         // Set Z boundary conditions
         if(eq.bcIds().bcs(eqId,fieldId).count(Dimensions::Simulation::SIM3D) > 0)
         {
            tau = Spectral::BoundaryConditions::tauMatrix(bound3D, eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM3D)->second);
            if(tau.first.nonZeros() > 0)
            {
               if(cZ != 1.0)
               {
                  tau.first *= cZ;
               }
               SparseMatrix tmp;
               Eigen::kroneckerProduct(tau.first, q1D, tmp);
               mat.first += tmp;
            }
            if(tau.second.nonZeros() > 0)
            {
               if(cZ != 1.0)
               {
                  tau.second *= cZ;
               }
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
