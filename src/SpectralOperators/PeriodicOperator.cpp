/** \file PeriodicOperator.cpp
 *  \brief Source of the implementation of the spectral operators with periodic directions
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
#include "SpectralOperators/PeriodicOperator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Spectral {

   MHDFloat PeriodicOperator::laplacian2D(const MHDFloat k, const MHDFloat m)
   {
      // Compute periodic laplacian
      MHDFloat lapl = (std::pow(k,2) + std::pow(m,2));

      return lapl;
   }

   MHDFloat PeriodicOperator::bilaplacian2D(const MHDFloat k, const MHDFloat m)
   {
      // Compute periodic bilaplacian
      MHDFloat lapl = (std::pow(k,4) + 2*std::pow(k,2)*std::pow(m,2) + std::pow(m,4));

      return lapl;
   }

   SparseMatrix PeriodicOperator::laplacian2D(const IOperator& op, const MHDFloat k, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = (op.diff(nBC, 2) - std::pow(k,2)*op.id(nBC));

      // Prune sparse laplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix PeriodicOperator::bilaplacian2D(const IOperator& op, const MHDFloat k, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = PeriodicOperator::laplacian2D(op, k, nBC);

      // Compute bilaplacian
      lapl = lapl*lapl;

      // Prune sparse bilaplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix PeriodicOperator::laplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = (op.diff(nBC, 2) - std::pow(k,2)*op.id(nBC) - std::pow(m,2)*op.id(nBC));

      // Prune sparse laplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix PeriodicOperator::laplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int nBCA, const int nBCB)
   {
      // Compute sparse laplacian
      SparseMatrix lapl;
      lapl = Eigen::kroneckerProduct(opB.diff(nBCB, 2), opA.id(nBCB));
      SparseMatrix tmp;
      tmp = Eigen::kroneckerProduct(opB.id(nBCB), opA.diff(nBCB,2));
      lapl = lapl + tmp;
      tmp = Eigen::kroneckerProduct(opB.id(nBCB), opA.id(nBCA));
      lapl = lapl - std::pow(k,2)*tmp;

      // Prune sparse laplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix PeriodicOperator::bilaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int nBC)
   {
      // Compute sparse laplacian
      SparseMatrix lapl = PeriodicOperator::laplacian3D(op, k, m, nBC);

      // Compute bilaplacian
      lapl = lapl*lapl;

      // Prune sparse bilaplacian
      lapl.prune(1e-32);

      return lapl;
   }

   SparseMatrix PeriodicOperator::qLaplacian2D(const IOperator& op, const MHDFloat k, const int p)
   {
      // Compute quasi inverse of laplacian
      SparseMatrix qLapl = (op.qDiff(p,2) - std::pow(k,2)*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

   SparseMatrix PeriodicOperator::qBilaplacian2D(const IOperator& op, const MHDFloat k, const int p)
   {
      // Compute quasi-inverse of bilaplacian
      SparseMatrix qLapl = (op.qDiff(p,4) - 2*std::pow(k,2)*op.qDiff(p,2) + std::pow(k,4)*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

   SparseMatrix PeriodicOperator::qLaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int pA, const int pB)
   {
      // Compute quasi inverse of laplacian
      SparseMatrix qLapl;
      qLapl = Eigen::kroneckerProduct(opB.qDiff(pB, 2), opA.qDiff(pA,0));
      SparseMatrix tmp;
      tmp = Eigen::kroneckerProduct(opB.qDiff(pB, 0), opA.qDiff(pA,2));
      qLapl = qLapl + tmp;
      tmp = Eigen::kroneckerProduct(opB.qDiff(pB, 0), opA.qDiff(pA,0));
      qLapl = qLapl - std::pow(k,2)*tmp;

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

   SparseMatrix PeriodicOperator::qBilaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int pA, const int pB)
   {
      // Compute quasi inverse of laplacian
      SparseMatrix qLapl;
      qLapl = Eigen::kroneckerProduct(opB.qDiff(pB,0), opA.qDiff(pA,4));
      SparseMatrix tmp;
      tmp = Eigen::kroneckerProduct(opB.qDiff(pB,4), opA.qDiff(pA,0));
      qLapl = qLapl + tmp;
      tmp = Eigen::kroneckerProduct(opB.qDiff(pB,2), opA.qDiff(pA,2));
      qLapl = qLapl + 2*std::pow(k,2)*tmp;
      tmp = Eigen::kroneckerProduct(opB.qDiff(pB,0), opA.qDiff(pA,2));
      qLapl = qLapl - 2*std::pow(k,2)*tmp;
      tmp = Eigen::kroneckerProduct(opB.qDiff(pB,2), opA.qDiff(pA,0));
      qLapl = qLapl - 2*std::pow(k,2)*tmp;
      tmp = Eigen::kroneckerProduct(opB.qDiff(pB,0), opA.qDiff(pA,0));
      qLapl = qLapl + std::pow(k,4)*tmp;

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

   SparseMatrix PeriodicOperator::qLaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int p)
   {
      // Compute quasi inverse of laplacian
      SparseMatrix qLapl = (op.qDiff(p,2) - (std::pow(k,2) - std::pow(m,2))*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

   SparseMatrix PeriodicOperator::qBilaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int p)
   {
      // Compute quasi-inverse of bilaplacian
      SparseMatrix qLapl = (op.qDiff(p,4) - (2*std::pow(k,2)+2*std::pow(m,2))*op.qDiff(p,2) + (std::pow(k,4) + 2*std::pow(k,2)*std::pow(m,2) + std::pow(m,4))*op.qDiff(p,0));

      // Prune sparse matrix
      qLapl.prune(1e-32);

      return qLapl;
   }

}
}
