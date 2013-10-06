/** 
 * @file EquationEigen1DTools.cpp
 * @brief Source of the tools for schemes with a single eigen direction
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <limits>

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Tools/EquationEigen1DTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

namespace Eigen1D {

   void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols)
   {
      // Get 1D dimension (fast)
      int nI = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM);
      // Get 2D dimension (slow)
      int nJ = spRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::TRANSFORM);

      // Get 3D dimension (medium)
      nMat = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();

      blocks.resize(nMat);
      blocks.setConstant(nI*nJ);
      cols.resize(nMat);
      cols.setConstant(1);
   }

   void computeKProduct(SparseMatrix& mat, const KRProduct& block)
   {
      mat = Eigen::kroneckerProduct(std::tr1::get<1>(block), std::tr1::get<0>(block));
   }

   void computeKProduct(DecoupledZSparse& mat, const KZProduct& block)
   {
      mat.real() = Eigen::kroneckerProduct(std::tr1::get<1>(block).real(), std::tr1::get<0>(block).real());
      mat.real() -= Eigen::kroneckerProduct(std::tr1::get<1>(block).imag(), std::tr1::get<0>(block).imag());
      mat.imag() = Eigen::kroneckerProduct(std::tr1::get<1>(block).real(), std::tr1::get<0>(block).imag());
      mat.imag() += Eigen::kroneckerProduct(std::tr1::get<1>(block).imag(), std::tr1::get<0>(block).real());
   }

   void computeKSum(SparseMatrix& mat, const KRSum& blocks)
   {
      if(blocks.size() > 0)
      {
         computeKProduct(mat, blocks.at(0));

         SparseMatrix tmp;
         for(std::vector<KRProduct>::const_iterator it = (blocks.begin()+1); it != blocks.end(); ++it)
         {
            computeKProduct(tmp, *it);

            mat += tmp;
         }
      }
   }

   void computeKSum(DecoupledZSparse& mat, const KZSum& blocks)
   {
      if(blocks.size() > 0)
      {
         computeKProduct(mat, blocks.at(0));

         DecoupledZSparse tmp;
         for(std::vector<KZProduct>::const_iterator it = (blocks.begin()+1); it != blocks.end(); ++it)
         {
            computeKProduct(tmp, *it);

            mat.real() += tmp.real();
            mat.imag() += tmp.imag();
         }
      }
   }

//   void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int p1D, const int p3D, const MHDFloat c1D, const MHDFloat c3D)
//   {
//      // Get 1D and 3D dimensions
//      int n1D = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
//      int n3D = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
//
//      // Create equation ID
//      SpectralFieldId eqId = std::make_pair(eq.name(), compId);
//
//      // Create spectral operators
//      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(n1D);
//      Spectral::OperatorSelector<Dimensions::Simulation::SIM3D>::Type spec3D(n3D);
//
//      // Create spectral boundary operators
//      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::BcType bound1D(n1D);
//      Spectral::OperatorSelector<Dimensions::Simulation::SIM3D>::BcType bound3D(n3D);
//
//      // Initialise output matrices
//      mat.real().resize(n1D*n3D,n1D*n3D);
//      mat.imag().resize(n1D*n3D,n1D*n3D);
//
//      // Nothing to do if c1D = 0 and c3D = 0;
//      if(c1D != 0 || c3D != 0)
//      {
//         // Storage for the boundary quasi-inverses
//         SparseMatrix q1D;
//         SparseMatrix q3D;
//
//         // Set boundary "operators"
//         if(eq.bcIds().hasEquation(eqId))
//         {
///// \mhdBug Number of boundary conditions might not be dealed with correctly!
//
//            // Impose 1D condition in corners
//            if(p1D == 0 && p3D > 0)
//            {
//               int nBC3D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM3D)->second.size();
//
//               // Set 1D boundary quasi-inverse
//               q1D = spec1D.id(p1D);
//
//               // Set 3D boundary quasi-inverse
//               q3D = spec3D.shiftId(nBC3D);
//
//            // Impose 3D condition in corners
//            } else if(p3D == 0 && p1D > 0)
//            {
//               int nBC1D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM1D)->second.size();
//
//               // Set 1D boundary quasi-inverse
//               q1D = spec1D.shiftId(nBC1D);
//
//               // Set 3D boundary quasi-inverse
//               q3D = spec3D.id(p3D);
//
//               // Special case for "no" boundary condition
//            } else if(p1D == 0 && p3D == 0 && (c1D == 0 || c3D == 0))
//            {
//               // Set 1D boundary quasi-inverse
//               q1D = spec1D.id(p1D);
//
//               // Set 3D boundary quasi-inverse
//               q3D = spec3D.id(p3D);
//            } else
//            {
//               throw Exception("Boundary conditions setup is incompatible!");
//            }
//            // Unknown equation
//         } else
//         {
//            throw Exception("Missing condition(s) for boundary operator!");
//         }
//         // Temporary storage for the boundary operators
//         DecoupledZSparse tau;
//
//         // Set boundary conditions on fieldId
//         if(eq.bcIds().hasField(eqId,fieldId))
//         {
//            // Set X boundary conditions
//            if(eq.bcIds().bcs(eqId,fieldId).count(Dimensions::Simulation::SIM1D) > 0)
//            {
//               tau = Spectral::BoundaryConditions::tauMatrix(bound1D, eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second);
//               if(tau.real().nonZeros() > 0)
//               {
//                  if(c1D != 1.0)
//                  {
//                     tau.real() *= c1D;
//                  }
//                  mat.real() = Eigen::kroneckerProduct(q3D, tau.real());
//               }
//
//               if(tau.imag().nonZeros() > 0)
//               {
//                  if(c1D != 1.0)
//                  {
//                     tau.imag() *= c1D;
//                  }
//                  mat.imag() = Eigen::kroneckerProduct(q3D, tau.imag());
//               }
//            }
//
//            // Set Z boundary conditions
//            if(eq.bcIds().bcs(eqId,fieldId).count(Dimensions::Simulation::SIM3D) > 0)
//            {
//               tau = Spectral::BoundaryConditions::tauMatrix(bound3D, eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM3D)->second);
//               if(tau.real().nonZeros() > 0)
//               {
//                  if(c3D != 1.0)
//                  {
//                     tau.real() *= c3D;
//                  }
//                  SparseMatrix tmp;
//                  tmp = Eigen::kroneckerProduct(tau.real(), q1D);
//                  mat.real() += tmp;
//               }
//               if(tau.imag().nonZeros() > 0)
//               {
//                  if(c3D != 1.0)
//                  {
//                     tau.imag() *= c3D;
//                  }
//                  SparseMatrix tmp;
//                  tmp = Eigen::kroneckerProduct(tau.imag(), q1D);
//                  mat.imag() += tmp;
//               }
//            }
//         }
//
//         // Prune matrices for safety
//         mat.real().prune(1e-32);
//         mat.imag().prune(1e-32);
//      }
//   }
}
}
}