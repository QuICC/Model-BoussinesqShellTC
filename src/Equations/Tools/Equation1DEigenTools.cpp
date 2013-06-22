/** \file Equation1DEigenTools.cpp
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
#include "Equations/Tools/Equation1DEigenTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   void Equation1DEigenTools::boundaryBlock1DEigen(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int p1D, const int p3D, const MHDFloat c1D, const MHDFloat c3D)
   {
      // Get 1D and 3D dimensions
      int n1D = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int n3D = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create equation ID
      SpectralFieldId eqId = std::make_pair(eq.name(), compId);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(n1D);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(n3D);

      // Create spectral boundary operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::BcType bound1D(n1D);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::BcType bound3D(n3D);

      // Initialise output matrices
      mat.first.resize(n1D*n3D,n1D*n3D);
      mat.second.resize(n1D*n3D,n1D*n3D);

      // Storage for the boundary quasi-inverses
      SparseMatrix q1D;
      SparseMatrix q3D;

      // Set boundary "operators"
      if(eq.bcIds().hasEquation(eqId))
      {
         // Set 1D boundary quasi-inverse
         q1D = spec1D.shiftId(p1D);

         // Set 3D boundary quasi-inverse
         q3D = spec3D.id(p3D);
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
               if(c1D != 1.0)
               {
                  tau.first *= c1D;
               }
               Eigen::kroneckerProduct(q3D, tau.first, mat.first);
            }

            if(tau.second.nonZeros() > 0)
            {
               if(c1D != 1.0)
               {
                  tau.second *= c1D;
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
               if(c3D != 1.0)
               {
                  tau.first *= c3D;
               }
               SparseMatrix tmp;
               Eigen::kroneckerProduct(tau.first, q1D, tmp);
               mat.first += tmp;
            }
            if(tau.second.nonZeros() > 0)
            {
               if(c3D != 1.0)
               {
                  tau.second *= c3D;
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
