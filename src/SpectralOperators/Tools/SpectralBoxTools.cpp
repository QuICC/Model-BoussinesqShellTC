/** 
 * @file SpectralBoxTools.cpp
 * @brief Source of the implementation of the spectral cartesian operators with (possibly) periodic directions
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "SpectralOperators/Tools/SpectralBoxTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   namespace Spectral {

      namespace BoxTools {

         KronEigen3DRSum laplacian2D(const MHDFloat k, const MHDFloat m)
         {
            KronEigen3DRSum lapl = -(std::pow(k,2) + std::pow(m,2));

            return lapl;
         }

         KronEigen3DRSum bilaplacian2D(const MHDFloat k, const MHDFloat m)
         {
            KronEigen3DRSum lapl = std::pow((std::pow(k,2) + std::pow(m,2)),2);

            return lapl;
         }

         KronEigen2DRSum laplacian2D(const IOperator& op, const MHDFloat k, const int nBC)
         {
            KronEigen2DRSum lapl = (op.diff(nBC, 2) - std::pow(k,2)*op.id(nBC));

            return lapl;
         }

         KronEigen2DRSum bilaplacian2D(const IOperator& op, const MHDFloat k, const int nBC)
         {
            KronEigen2DRSum lapl = (op.diff(nBC, 4) - 2.0*std::pow(k,2)*op.diff(nBC,2) + std::pow(k,4)*op.id(nBC));

            return lapl;
         }

         KronEigen1DRSum laplacian2D(const IOperator& opA, const IOperator& opB, const int nBCA, const int nBCB)
         {
            KronEigen1DRSum lapl;
            lapl.push_back(std::tr1::make_tuple(opA.id(nBCA), opB.diff(nBCB, 2)));
            lapl.push_back(std::tr1::make_tuple(opA.diff(nBCA, 2), opB.id(nBCB)));

            return lapl;
         }

         KronEigen1DRSum bilaplacian2D(const IOperator& opA, const IOperator& opB, const int nBCA, const int nBCB)
         {
            KronEigen1DRSum lapl;
            lapl.push_back(std::tr1::make_tuple(opA.id(nBCA), opB.diff(nBCB, 4)));
            lapl.push_back(std::tr1::make_tuple(2.0*opA.diff(nBCA,2), opB.diff(nBCB, 2)));
            lapl.push_back(std::tr1::make_tuple(opA.diff(nBCA, 4), opB.id(nBCB)));

            return lapl;
         }

         KronEigen3DRSum laplacian3D(const MHDFloat k, const MHDFloat l, const MHDFloat m)
         {
            KronEigen3DRSum lapl = -(std::pow(k,2) + std::pow(l,2) + std::pow(m,2));

            return lapl;
         }

         KronEigen3DRSum bilaplacian3D(const MHDFloat k, const MHDFloat l, const MHDFloat m)
         {
            KronEigen3DRSum lapl = std::pow((std::pow(k,2) + std::pow(l,2) + std::pow(m,2)),2);

            return lapl;
         }

         KronEigen2DRSum laplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int nBC)
         {
            KronEigen2DRSum lapl = (op.diff(nBC, 2) - (std::pow(k,2) + std::pow(m,2))*op.id(nBC));

            return lapl;
         }

         KronEigen2DRSum bilaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int nBC)
         {
            MHDFloat c = std::pow(k,2) + std::pow(m,2);
            KronEigen2DRSum lapl = (op.diff(nBC, 4) - 2.0*c*op.diff(nBC,2) + std::pow(c,2)*op.id(nBC));

            return lapl;
         }

         KronEigen1DRSum laplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int nBCA, const int nBCB)
         {
            KronEigen1DRSum lapl;
            lapl.push_back(std::tr1::make_tuple(opA.id(nBCA), opB.diff(nBCB, 2)));
            lapl.push_back(std::tr1::make_tuple(-std::pow(k,2)*opA.id(nBCA), opB.id(nBCB)));
            lapl.push_back(std::tr1::make_tuple(opA.diff(nBCA,2), opB.id(nBCB)));

            return lapl;
         }

         KronEigen1DRSum bilaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int nBCA, const int nBCB)
         {
            KronEigen1DRSum lapl;
            lapl.push_back(std::tr1::make_tuple(opA.id(nBCA), opB.diff(nBCB, 4)));
            lapl.push_back(std::tr1::make_tuple(-2.0*std::pow(k,2)*opA.id(nBCA), opB.diff(nBCB,2)));
            lapl.push_back(std::tr1::make_tuple(2.0*opA.diff(nBCA,2), opB.diff(nBCB, 2)));
            lapl.push_back(std::tr1::make_tuple(std::pow(k,4)*opA.id(nBCA), opB.id(nBCB)));
            lapl.push_back(std::tr1::make_tuple(-2.0*std::pow(k,2)*opA.diff(nBCA,2), opB.id(nBCB)));
            lapl.push_back(std::tr1::make_tuple(opA.diff(nBCA,4), opB.id(nBCB)));

            return lapl;
         }

         KronNoEigenRSum laplacian3D(const IOperator& opA, const IOperator& opB, const IOperator& opC, const int nBCA, const int nBCB, const int nBCC)
         {
            KronNoEigenRSum lapl;
            lapl.push_back(std::tr1::make_tuple(opA.diff(nBCA,2), opB.id(nBCB), opC.id(nBCC)));
            lapl.push_back(std::tr1::make_tuple(opA.id(nBCA), opB.diff(nBCB,2), opC.id(nBCC)));
            lapl.push_back(std::tr1::make_tuple(opA.id(nBCA), opB.id(nBCB), opC.diff(nBCC,2)));

            return lapl;
         }

         KronNoEigenRSum bilaplacian3D(const IOperator& opA, const IOperator& opB, const IOperator& opC, const int nBCA, const int nBCB, const int nBCC)
         {
            KronNoEigenRSum lapl;
            lapl.push_back(std::tr1::make_tuple(2.0*opA.diff(nBCA,2), opB.diff(nBCB,2), opC.id(nBCC)));
            lapl.push_back(std::tr1::make_tuple(2.0*opA.id(nBCA), opB.diff(nBCB,2), opC.diff(nBCC,2)));
            lapl.push_back(std::tr1::make_tuple(2.0*opA.diff(nBCA,2), opB.id(nBCB), opC.diff(nBCC,2)));
            lapl.push_back(std::tr1::make_tuple(opA.diff(nBCA,4), opB.id(nBCB), opC.id(nBCC)));
            lapl.push_back(std::tr1::make_tuple(opA.id(nBCA), opB.diff(nBCB,4), opC.id(nBCC)));
            lapl.push_back(std::tr1::make_tuple(opA.id(nBCA), opB.id(nBCB), opC.diff(nBCC,4)));

            return lapl;
         }

         KronEigen2DRSum qLaplacian2D(const IOperator& op, const MHDFloat k, const int p)
         {
            KronEigen2DRSum qLapl = (op.qDiff(p,2) - std::pow(k,2)*op.qDiff(p,0));

            return qLapl;
         }

         KronEigen2DRSum qBilaplacian2D(const IOperator& op, const MHDFloat k, const int p)
         {
            KronEigen2DRSum qLapl = (op.qDiff(p,4) - 2*std::pow(k,2)*op.qDiff(p,2) + std::pow(k,4)*op.qDiff(p,0));

            return qLapl;
         }

         KronEigen1DRSum qLaplacian2D(const IOperator& opA, const IOperator& opB, const int pA, const int pB)
         {
            KronEigen1DRSum qLapl;
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA, 0), opB.qDiff(pB, 2)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA, 2), opB.qDiff(pB,0)));

            return qLapl;
         }

         KronEigen1DRSum qBilaplacian2D(const IOperator& opA, const IOperator& opB, const int pA, const int pB)
         {
            KronEigen1DRSum qLapl;
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,0), opB.qDiff(pB, 4)));
            qLapl.push_back(std::tr1::make_tuple(2.0*opA.qDiff(pA,2), opB.qDiff(pB, 2)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA, 4), opB.qDiff(pB,0)));

            return qLapl;
         }

         KronEigen2DRSum qLaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int p)
         {
            SparseMatrix qLapl = (op.qDiff(p,2) - (std::pow(k,2) + std::pow(m,2))*op.qDiff(p,0));

            return qLapl;
         }

         KronEigen2DRSum qBilaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int p)
         {
            MHDFloat c = std::pow(k,2) + std::pow(m,2);
            KronEigen2DRSum qLapl = (op.qDiff(p, 4) - 2.0*c*op.qDiff(p,2) + std::pow(c,2)*op.qDiff(p,0));

            return qLapl;
         }

         KronEigen1DRSum qLaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int pA, const int pB)
         {
            KronEigen1DRSum qLapl;
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,0), opB.qDiff(pB, 2)));
            qLapl.push_back(std::tr1::make_tuple(-std::pow(k,2)*opA.qDiff(pA, 0), opB.qDiff(pB,0)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,2), opB.qDiff(pB,0)));

            return qLapl;
         }

         KronEigen1DRSum qBilaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int pA, const int pB)
         {
            KronEigen1DRSum qLapl;
            qLapl.push_back(std::tr1::make_tuple(-2.0*std::pow(k,2)*opA.qDiff(pA,0), opB.qDiff(pB,2)));
            qLapl.push_back(std::tr1::make_tuple(2.0*opA.qDiff(pA,2), opB.qDiff(pB, 2)));
            qLapl.push_back(std::tr1::make_tuple(-2.0*std::pow(k,2)*opA.qDiff(pA,2), opB.qDiff(pB,0)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,0), opB.qDiff(pB, 4)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,4), opB.qDiff(pB,0)));
            qLapl.push_back(std::tr1::make_tuple(std::pow(k,4)*opA.qDiff(pA,0), opB.qDiff(pB,0)));

            return qLapl;
         }

         KronNoEigenRSum qLaplacian3D(const IOperator& opA, const IOperator& opB, const IOperator& opC, const int pA, const int pB, const int pC)
         {
            KronNoEigenRSum qLapl;
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,2), opB.qDiff(pB,0), opC.qDiff(pC,0)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,0), opB.qDiff(pB,2), opC.qDiff(pC,0)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,0), opB.qDiff(pB,0), opC.qDiff(pC,2)));

            return qLapl;
         }

         KronNoEigenRSum qBilaplacian3D(const IOperator& opA, const IOperator& opB, const IOperator& opC, const int pA, const int pB, const int pC)
         {
            KronNoEigenRSum qLapl;
            qLapl.push_back(std::tr1::make_tuple(2.0*opA.qDiff(pA,2), opB.qDiff(pB,2), opC.qDiff(pC,0)));
            qLapl.push_back(std::tr1::make_tuple(2.0*opA.qDiff(pA,0), opB.qDiff(pB,2), opC.qDiff(pC,2)));
            qLapl.push_back(std::tr1::make_tuple(2.0*opA.qDiff(pA,2), opB.qDiff(pB,0), opC.qDiff(pC,2)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,4), opB.qDiff(pB,0), opC.qDiff(pC,0)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,0), opB.qDiff(pB,4), opC.qDiff(pC,0)));
            qLapl.push_back(std::tr1::make_tuple(opA.qDiff(pA,0), opB.qDiff(pB,0), opC.qDiff(pC,4)));

            return qLapl;
         }

      }
   }
}
