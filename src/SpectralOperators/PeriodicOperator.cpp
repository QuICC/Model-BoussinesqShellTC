/** \file PeriodicOperator.cpp
 *  \brief Source of the implementation of the spectral operator for the Chebyshev basis
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "SpectralOperators/PeriodicOperator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Spectral {

   SparseMatrix PeriodicOperator::laplacian2D(const IOperator& op, const int k, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = (op.diff(nBC, 2) - static_cast<MHDFloat>(std::pow(k,2))*op.id(nBC));

      // Prune sparse laplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix PeriodicOperator::bilaplacian2D(const IOperator& op, const int k, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = PeriodicOperator::laplacian2D(op, k, nBC);

      // Compute bilaplacian
      lapl = lapl*lapl;

      // Prune sparse bilaplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix PeriodicOperator::laplacian3D(const IOperator& op, const int k, const int m, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = (op.diff(nBC, 2) - static_cast<MHDFloat>(std::pow(k,2))*op.id(nBC) - static_cast<MHDFloat>(std::pow(m,2))*op.id(nBC));

      // Prune sparse laplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix PeriodicOperator::bilaplacian3D(const IOperator& op, const int k, const int m, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = PeriodicOperator::laplacian3D(op, k, m, nBC);

      // Compute bilaplacian
      lapl = lapl*lapl;

      // Prune sparse bilaplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix PeriodicOperator::qLaplacian2D(const IOperator& op, const int k, const int p)
   {
      // Ccompute quasi inverse of laplacian
      SparseMatrix qLapl = (op.qDiff(p,2) - static_cast<MHDFloat>(std::pow(k,2))*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

   SparseMatrix PeriodicOperator::qBilaplacian2D(const IOperator& op, const int k, const int p)
   {
      // Compute quasi-inverse of bilaplacian
      SparseMatrix qLapl = (op.qDiff(p,4) - static_cast<MHDFloat>(2*std::pow(k,2))*op.qDiff(p,2) + static_cast<MHDFloat>(std::pow(k,4))*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

   SparseMatrix PeriodicOperator::qLaplacian3D(const IOperator& op, const int k, const int m, const int p)
   {
      // Ccompute quasi inverse of laplacian
      SparseMatrix qLapl = (op.qDiff(p,2) - static_cast<MHDFloat>(std::pow(k,2))*op.qDiff(p,0) - static_cast<MHDFloat>(std::pow(m,2))*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

   SparseMatrix PeriodicOperator::qBilaplacian3D(const IOperator& op, const int k, const int m, const int p)
   {
      // Compute quasi-inverse of bilaplacian
      SparseMatrix qLapl = (op.qDiff(p,4) - static_cast<MHDFloat>(2*std::pow(k,2)+2*std::pow(m,2))*op.qDiff(p,2) + static_cast<MHDFloat>(std::pow(k,4) + 2*std::pow(k,2)*std::pow(m,2) + std::pow(m,4))*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

}
}
