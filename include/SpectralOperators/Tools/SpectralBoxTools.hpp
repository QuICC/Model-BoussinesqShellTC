/** 
 * @file SpectralBoxTools.hpp
 * @brief Implementation of the multidimensional cartesian box operators with possibly periodic dimensions
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPECTRALBOXTOOLS_HPP
#define SPECTRALBOXTOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/KroneckerTypedefs.hpp"
#include "SpectralOperators/IOperator.hpp"

namespace GeoMHDiSCC {

   namespace Spectral {

      namespace BoxTools {

         /**
          * @brief Compute the 2D Laplacian operator for 2D periodic
          *
          * @param k    First wave number
          * @param m    Second wave number
          */
         KronEigen3DRSum laplacian2D(const MHDFloat k, const MHDFloat m);

         /**
          * @brief Compute the 2D bi-Laplacian operator for 2D periodic
          *
          * @param k    First wave number
          * @param m    Second wave number
          */
         KronEigen3DRSum bilaplacian2D(const MHDFloat k, const MHDFloat m);

         /**
          * @brief Compute the 2D Laplacian operator for 1D periodic
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param nBC  Number of boundary rows
          */
         KronEigen2DRSum laplacian2D(const IOperator& op, const MHDFloat k, const int nBC);

         /**
          * @brief Compute the 2D bi-Laplacian operator for 1D periodic
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param nBC  Number of boundary rows
          */
         KronEigen2DRSum bilaplacian2D(const IOperator& op, const MHDFloat k, const int nBC);

         /**
          * @brief Compute the 2D Laplacian operator
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param nBCA  Number of boundary rows for first operator
          * @param nBCB  Number of boundary rows for second operator
          */
         KronEigen1DRSum laplacian2D(const IOperator& opA, const IOperator& opB, const int nBCA, const int nBCB);

         /**
          * @brief Compute the 2D bi-Laplacian operator
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param nBCA Number of boundary rows for first operator
          * @param nBCB x Number of boundary rows for second operator
          */
         KronEigen1DRSum bilaplacian2D(const IOperator& opA, const IOperator& opB, const int nBCA, const int nBCB);

         /**
          * @brief Compute the 3D Laplacian operator for 3D periodic
          *
          * @param k First wave number
          * @param l Second wave number
          * @param m Third wave number operator
          */
         KronEigen3DRSum laplacian3D(const MHDFloat k, const MHDFloat l, const MHDFloat m);

         /**
          * @brief Compute the 3D bi-Laplacian operator for 3D periodic
          *
          * @param k First wave number
          * @param l Second wave number
          * @param m Third wave number operator
          */
         KronEigen3DRSum bilaplacian3D(const MHDFloat k, const MHDFloat l, const MHDFloat m);

         /**
          * @brief Compute the 3D Laplacian operator for 2D periodic
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param m    Second wave number
          * @param nBC  Number of boundary rows
          */
         KronEigen2DRSum laplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int nBC);

         /**
          * @brief Compute the 3D bi-Laplacian operator for 2D periodic
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param m    Second wave number
          * @param nBC  Number of boundary rows
          */
         KronEigen2DRSum bilaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int nBC);

         /**
          * @brief Compute the 3D Laplacian operator for 1D periodic
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param k    First wave number
          * @param nBCA Number of boundary rows for first operator
          * @param nBCB Number of boundary rows for second operator
          */
         KronEigen1DRSum laplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int nBCA, const int nBCB);

         /**
          * @brief Compute the 3D bi-Laplacian operator for 1D periodic
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param k    First wave number
          * @param nBCA Number of boundary rows for first operator
          * @param nBCB Number of boundary rows for second operator
          */
         KronEigen1DRSum bilaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int nBCA, const int nBCB);

         /**
          * @brief Compute the 3D Laplacian operator
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param opC  Third spectral operator
          * @param nBCA Number of boundary rows for first operator
          * @param nBCB Number of boundary rows for second operator
          * @param nBCC Number of boundary rows for third operator
          */
         KronNoEigenRSum laplacian3D(const IOperator& opA, const IOperator& opB, const IOperator& opC, const int nBCA, const int nBCB, const int nBCC);

         /**
          * @brief Compute the 3D bi-Laplacian operator
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param opC  Third spectral operator
          * @param nBCA Number of boundary rows for first operator
          * @param nBCB Number of boundary rows for second operator
          * @param nBCC Number of boundary rows for third operator
          */
         KronNoEigenRSum bilaplacian3D(const IOperator& opA, const IOperator& opB, const IOperator& opC, const int nBCA, const int nBCB, const int nBCC);

         /**
          * @brief Compute the 2D quasi inverse Laplacian operator for 1D periodic
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param p    Order of the quasi inverse
          */
         KronEigen2DRSum qLaplacian2D(const IOperator& op, const MHDFloat k, const int p);

         /**
          * @brief Compute the 2D quasi inverse bi-Laplacian operator  for 1D periodic
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param p    Order of the quasi inverse
          */
         KronEigen2DRSum qBilaplacian2D(const IOperator& op, const MHDFloat k, const int p);

         /**
          * @brief Compute the 2D quasi inverse Laplacian operator
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param pA   Order of the quasi inverse for first operator
          * @param pB   Order of the quasi inverse for second operator
          */
         KronEigen1DRSum qLaplacian2D(const IOperator& opA, const IOperator& opB, const int pA, const int pB);

         /**
          * @brief Compute the 2D quasi inverse bi-Laplacian operator
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param pA   Order of the quasi inverse for first operator
          * @param pB   Order of the quasi inverse for second operator
          */
         KronEigen1DRSum qBilaplacian2D(const IOperator& opA, const IOperator& opB, const int pA, const int pB);

         /**
          * @brief Compute the 3D quasi inverse Laplacian operator in 2D periodic
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param m    Second wave number
          * @param p    Order of the quasi inverse
          */
         KronEigen2DRSum qLaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int p);

         /**
          * @brief Compute the 3D quasi inverse bi-Laplacian operator in 2D periodic
          *
          * @param op   Spectral operator
          * @param k    First wave number
          * @param m    Second wave number
          * @param p    Order of the quasi inverse
          */
         KronEigen2DRSum qBilaplacian3D(const IOperator& op, const MHDFloat k, const MHDFloat m, const int p);

         /**
          * @brief Compute the 3D quasi inverse Laplacian operator in 1D periodic
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param k    First wave number
          * @param pA   Order of the quasi inverse for first operator
          * @param pB   Order of the quasi inverse for second operator
          */
         KronEigen1DRSum qLaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int pA, const int pB);

         /**
          * @brief Compute the 3D quasi inverse bi-Laplacian operator in 1D periodic 
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param k    First wave number
          * @param pA   Order of the quasi inverse for first operator
          * @param pB   Order of the quasi inverse for second operator
          */
         KronEigen1DRSum qBilaplacian3D(const IOperator& opA, const IOperator& opB, const MHDFloat k, const int pA, const int pB);

         /**
          * @brief Compute the 3D quasi inverse Laplacian operator
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param opC  Third spectral operator
          * @param pA   Order of the quasi inverse for first operator
          * @param pB   Order of the quasi inverse for second operator
          * @param pC   Order of the quasi inverse for third operator
          */
         KronNoEigenRSum qLaplacian3D(const IOperator& opA, const IOperator& opB, const IOperator& opC, const int pA, const int pB, const int pC);

         /**
          * @brief Compute the 3D quasi inverse bi-Laplacian operator
          *
          * @param opA  First spectral operator
          * @param opB  Second spectral operator
          * @param opC  Third spectral operator
          * @param pA   Order of the quasi inverse for first operator
          * @param pB   Order of the quasi inverse for second operator
          * @param pC   Order of the quasi inverse for third operator
          */
         KronNoEigenRSum qBilaplacian3D(const IOperator& opA, const IOperator& opB, const IOperator& opC, const int pA, const int pB, const int pC);
      }
   }
}

#endif // SPECTRALBOXTOOLS_HPP
