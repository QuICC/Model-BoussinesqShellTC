/** 
 * @file Equation3DEigenTools.cpp
 * @brief Source of the tools for schemes with three eigen direction
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
#include "Equations/Tools/Equation3DEigenTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   const Boundary::BCIndex Equation3DEigenTools::INDEPENDENT = std::tr1::make_tuple(std::numeric_limits<int>::min(),std::numeric_limits<int>::min(),std::numeric_limits<int>::min());

   void Equation3DEigenTools::makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols)
   {
      // Get 1D dimension (fast)
      int nI = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM);

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i)*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATB1D>(i,j);
      }

      blocks.resize(nMat);
      blocks.setConstant(1);
      cols.resize(nMat);
      cols.setConstant(1);
   }

   template <typename TEquation> void Equation3DEigenTools::storeBoundaryCondition(IEquation& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& coeffs, const std::vector<Equation3DEigenTools::BCIndex>& bcIdx)
   {
      assert(coeffs.size() == isEigsDependent.size());
      assert(coeffs.size() == 0);
   }

//   void Equation3DEigenTools::boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat c1D)
//   {
//      throw Exception("Not yet implemented!");
//   }
}
}
