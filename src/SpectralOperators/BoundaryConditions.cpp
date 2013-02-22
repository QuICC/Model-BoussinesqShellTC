/** \file BoundaryConditions.cpp
 *  \brief Source of the implementation of the spectral boundary conditions
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "SpectralOperators/BoundaryConditions.hpp"

// Project includes
//
#include "Exception/Exception.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   BoundaryConditions::BoundaryConditions()
   {
   }

   BoundaryConditions::~BoundaryConditions()
   {
   }

   DecoupledZMatrix BoundaryConditions::tauLines(const IBoundary& bcOp, const std::map<BoundaryConditions::Id,IBoundary::Position>& bcId)
   {
      // Count boundary conditions
      int nBCs = bcId.size();

      // Flags to know if real and imaginary parts are used
      bool hasReal = false;
      bool hasImag = false;

      // Check for compatible sizes
      assert(bcOp.basisN() >= nBCs);

      // Initialise tau lines matrix
      DecoupledZMatrix tauLines(std::make_pair(Matrix(bcOp.basisN(), nBCs),Matrix(bcOp.basisN(), nBCs)));
      tauLines.first.setZero();
      tauLines.second.setZero();

      // Map iterator
      std::map<BoundaryConditions::Id,BoundaryConditions::Position>::const_iterator mapIt;

      // Create boundary values
      int idx = 0;
      for(mapIt = bcId.begin(); mapIt != bcId.end(); ++mapIt)
      {
         switch(mapIt->first)
         {
            case BoundaryConditions::VALUE:
               tauLines.first.col(idx) = bcOp.value(mapIt->second);
               idx++;
               hasReal = true;
               break;
            case BoundaryConditions::FIRST_DERIVATIVE:
               tauLines.first.col(idx) = bcOp.firstDerivative(mapIt->second);
               idx++;
               hasReal = true;
               break;
            case BoundaryConditions::SECOND_DERIVATIVE:
               tauLines.first.col(idx) = bcOp.secondDerivative(mapIt->second);
               idx++;
               hasReal = true;
               break;
            case BoundaryConditions::BETA_SLOPE:
               tauLines.second.col(idx) = static_cast<MHDFloat>(-k)*bcOp.value(mapIt->second);
               idx++;
               hasImag = true;
               break;
            default:
               throw Exception("Unknown boundary condition for Chebyshev tau method");
               break;
         }
      }

      // Clear real matrix if not conditions have been set
      if(!hasReal)
      {
         tauLines.first.resize(0,0);
      }

      // Clear imaginary matrix if not conditions have been set
      if(!hasImag)
      {
         tauLines.second.resize(0,0);
      }

      return tauLines;
   }

}
}
