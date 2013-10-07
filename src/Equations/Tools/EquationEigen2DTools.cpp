/** 
 * @file EquationEigen2DTools.cpp
 * @brief Source of the tools for schemes with two eigen direction
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
#include "Equations/Tools/EquationEigen2DTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

namespace Eigen2D {

   void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols)
   {
      // Get 1D dimension (fast)
      int nI = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM);

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      blocks.resize(nMat);
      blocks.setConstant(nI);
      cols.resize(nMat);
      cols.setConstant(1);
   }

   void computeKProduct(SparseMatrix& mat, const KRProduct& block)
   {
      assert(block.size() > 0);

      mat = block;
   }

   void computeKProduct(DecoupledZSparse& mat, const KZProduct& block)
   {
      assert(block.real().size() > 0);
      assert(block.real().size() == block.imag().size());

      mat.real() = block.real();
      mat.imag() = block.imag();
   }

   void computeKSum(SparseMatrix& mat, const KRSum& blocks)
   {
      computeKProduct(mat, blocks);
   }

   void computeKSum(DecoupledZSparse& mat, const KZSum& blocks)
   {
      computeKProduct(mat, blocks);
   }

//   void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat c1D)
//   {
//      // Get 1D dimension
//      int n1D = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
//
//      // Create equation ID
//      SpectralFieldId eqId = std::make_pair(eq.name(), compId);
//
//      // Create spectral operators
//      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(n1D);
//
//      // Create spectral boundary operators
//      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::BcType bound1D(n1D);
//
//      // Initialise output matrices
//      mat.real().resize(n1D,n1D);
//      mat.imag().resize(n1D,n1D);
//
//      if(c1D != 0)
//      {
//         // Safety check to abort early
//         if(!eq.bcIds().hasEquation(eqId))
//         {
//            throw Exception("Missing condition(s) for boundary operator!");
//         }
//
//         // Set boundary conditions on fieldId
//         if(eq.bcIds().hasField(eqId,fieldId))
//         {
//            // Set 1D boundary conditions
//            if(eq.bcIds().bcs(eqId,fieldId).count(Dimensions::Simulation::SIM1D) > 0)
//            {
//               mat = Spectral::BoundaryConditions::tauMatrix(bound1D, eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second);
//               if(mat.real().nonZeros() > 0)
//               {
//                  if(c1D != 1.0)
//                  {
//                     mat.real() *= c1D;
//                  }
//               }
//
//               if(mat.imag().nonZeros() > 0)
//               {
//                  if(c1D != 1.0)
//                  {
//                     mat.imag() *= c1D;
//                  }
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
