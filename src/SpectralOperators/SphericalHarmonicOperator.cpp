/** 
 * @file SphericalHarmonicOperator.cpp
 * @brief Source of the implementation of the spherical harmonic spectral operators
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// System includes
//
#include <cassert>

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "SpectralOperators/SphericalHarmonicOperator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Spectral {

   MHDFloat SphericalHarmonicOperator::surfaceLaplacian(const MHDFloat l, const MHDFloat m)
   {
      // Compute periodic laplacian
      MHDFloat lapl = -l*(l+1);

      return lapl;
   }

   MHDFloat SphericalHarmonicOperator::surfaceBilaplacian(const MHDFloat l, const MHDFloat m)
   {
      // Compute periodic bilaplacian
      MHDFloat lapl = l*(l+1)*l*(l+1);

      return lapl;
   }

   SparseMatrix SphericalHarmonicOperator::laplacian(const IOperator& op, const MHDFloat l, const MHDFloat m, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = (op.diff(nBC, 2) + SphericalHarmonicOperator::surfaceLaplacian(l,m)*op.id(nBC));

      // Prune sparse laplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix SphericalHarmonicOperator::bilaplacian(const IOperator& op, const MHDFloat l, const MHDFloat m, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = (op.diff(nBC, 4) + 2*SphericalHarmonicOperator::surfaceLaplacian(l,m)*op.diff(nBC,2) + SphericalHarmonicOperator::surfaceBilaplacian(l,m)*op.id(nBC));

      // Prune sparse bilaplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix SphericalHarmonicOperator::qLaplacian(const IOperator& op, const MHDFloat l, const MHDFloat m, const int p)
   {
      // Compute quasi inverse of laplacian
      SparseMatrix qLapl = (op.qDiff(p,2) + SphericalHarmonicOperator::surfaceLaplacian(l,m)*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

   SparseMatrix SphericalHarmonicOperator::qBilaplacian(const IOperator& op, const MHDFloat l, const MHDFloat m, const int p)
   {
      // Compute quasi-inverse of bilaplacian
      SparseMatrix qLapl = (op.qDiff(p,4) + 2*SphericalHarmonicOperator::surfaceLaplacian(l,m)*op.qDiff(p,2) + SphericalHarmonicOperator::surfaceBilaplacian(l,m)*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

}
}
