/** 
 * @file Equation1DEigenTools.cpp
 * @brief Source of the tools for schemes with a single eigen direction
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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

   void Equation1DEigenTools::boundaryBlock1DEigen(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int p1D, const int p3D, const MHDFloat c1D, const MHDFloat c3D)
   {
      // Get 1D and 3D dimensions
      int n1D = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int n3D = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

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

      // Nothing to do if c1D = 0 and c3D = 0;
      if(c1D != 0 || c3D != 0)
      {
         // Storage for the boundary quasi-inverses
         SparseMatrix q1D;
         SparseMatrix q3D;

         // Set boundary "operators"
         if(eq.bcIds().hasEquation(eqId))
         {
/// \mhdBug Number of boundary conditions might not be dealed with correctly!

            // Impose 1D condition in corners
            if(p1D == 0 && p3D > 0)
            {
               int nBC3D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM3D)->second.size();

               // Set 1D boundary quasi-inverse
               q1D = spec1D.id(p1D);

               // Set 3D boundary quasi-inverse
               q3D = spec3D.shiftId(nBC3D);

            // Impose 3D condition in corners
            } else if(p3D == 0 && p1D > 0)
            {
               int nBC1D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM1D)->second.size();

               // Set 1D boundary quasi-inverse
               q1D = spec1D.shiftId(nBC1D);

               // Set 3D boundary quasi-inverse
               q3D = spec3D.id(p3D);

               // Special case for "no" boundary condition
            } else if(p1D == 0 && p3D == 0 && (c1D == 0 || c3D == 0))
            {
               // Set 1D boundary quasi-inverse
               q1D = spec1D.id(p1D);

               // Set 3D boundary quasi-inverse
               q3D = spec3D.id(p3D);
            } else
            {
               throw Exception("Boundary conditions setup is incompatible!");
            }
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
                  mat.first = Eigen::kroneckerProduct(q3D, tau.first);
               }

               if(tau.second.nonZeros() > 0)
               {
                  if(c1D != 1.0)
                  {
                     tau.second *= c1D;
                  }
                  mat.second = Eigen::kroneckerProduct(q3D, tau.second);
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
                  tmp = Eigen::kroneckerProduct(tau.first, q1D);
                  mat.first += tmp;
               }
               if(tau.second.nonZeros() > 0)
               {
                  if(c3D != 1.0)
                  {
                     tau.second *= c3D;
                  }
                  SparseMatrix tmp;
                  tmp = Eigen::kroneckerProduct(tau.second, q1D);
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
}
