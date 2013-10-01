/** 
 * @file EquationNoEigenTools.cpp
 * @brief Source of the tools for schemes with no eigen direction
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Equations/Tools/EquationNoEigenTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   void EquationNoEigenTools::makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols)
   {
      // Get 1D dimension (fast)
      int nI = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM);
      // Get 2D dimension (medium)
      int nK = spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::TRANSFORM);
      // Get 3D dimension (slow)
      int nJ = spRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::TRANSFORM);

      // Get 3D dimension (medium)
      nMat = 1;

      blocks.resize(nMat);
      blocks.setConstant(nI*nJ*nK);
      cols.resize(nMat);
      cols.setConstant(1);
   }

   void EquationNoEigenTools::boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int p1D, const int p3D, const MHDFloat c1D, const MHDFloat c3D)
   {
      throw Exception("Not yet implemented!");
   }
}
}
