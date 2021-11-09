/**
 * @file WorlandOperators.cpp
 * @brief Source of the implementation of the Jones-Worland operators
 */

// System includes
//

// External includes
//

// Class include
//
#include "PolynomialTransforms/WorlandOperators.hpp"

// Project includes
//
#include "Base/Precision.hpp"
#include "Exceptions/Exception.hpp"
#include "Quadratures/WorlandChebyshevRule.hpp"
#include "Quadratures/LegendreRule.hpp"
#include "PolynomialTransforms/WorlandPolynomial.hpp"

namespace QuICC {

namespace Polynomial {

   void WorlandOperators::integralR3(Matrix& op, const int nN, const int l)
   {
      Matrix tOp, tOp2;

      // Compute integral weights
      int specSize = nN + 2;
      int gridSize = std::lround(std::ceil(4.0*nN/3.) + 1);
      Array legGrid(gridSize);
      Array legWeights(legGrid.size());
      internal::Array ilegGrid, ilegWeights;
      LegendreRule::computeQuadrature(legGrid, legWeights, ilegGrid, ilegWeights, legGrid.size());
      legGrid.array() = ((legGrid.array() + 1.0)/2.0);
      ilegGrid.array() = ((ilegGrid.array() + 1.0)/2.0);

      tOp.resize(legGrid.size(), specSize);
      internal::Matrix  itmp, ipoly;
      Polynomial::WorlandPolynomial::Wnl(tOp, ipoly, l+1, ilegGrid);
      Array intgWeights = 0.5*(legWeights.asDiagonal()*tOp).colwise().sum().transpose();

      Array chebGrid(gridSize);
      Array chebWeights(chebGrid.size());
      internal::Array ichebGrid, ichebWeights;
      WorlandChebyshevRule::computeQuadrature(chebGrid, chebWeights, ichebGrid, ichebWeights, gridSize);

      tOp.resize(legGrid.size(), nN);
      tOp2.resize(legGrid.size(), specSize);
      Polynomial::WorlandPolynomial::Wnl(tOp, itmp, l, ichebGrid);
      Polynomial::WorlandPolynomial::Wnl(tOp2, itmp, l+1, ichebGrid);
      Array r3 = chebWeights.array()*chebGrid.array().pow(3);
      op = (intgWeights.topRows(specSize).transpose()*tOp2.transpose()*r3.asDiagonal()*tOp).transpose();
   }

}
}
