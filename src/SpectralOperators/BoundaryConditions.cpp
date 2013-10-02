/** 
 * @file BoundaryConditions.cpp
 * @brief Source of the implementation of the spectral boundary conditions
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "SpectralOperators/BoundaryConditions.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   BoundaryConditions::BoundaryConditions()
   {
   }

   BoundaryConditions::~BoundaryConditions()
   {
   }

   DecoupledZMatrix BoundaryConditions::tauLines(const IBoundary& bcOp, const std::vector<std::pair<BoundaryConditions::Id,IBoundary::Position> >& bcId)
   {
      // Count boundary conditions
      int nBCs = bcId.size();

      // Flags to know if real and imaginary parts are used
      bool hasReal = false;
      bool hasImag = false;

      // Check for compatible sizes
      assert(bcOp.basisN() >= nBCs);

      // Initialise tau lines matrix
      DecoupledZMatrix lines(bcOp.basisN(), nBCs);
      lines.setZero();

      // Map iterator
      std::vector<std::pair<BoundaryConditions::Id,IBoundary::Position> >::const_iterator it;

      // Create boundary values
      int idx = 0;
      for(it = bcId.begin(); it != bcId.end(); ++it)
      {
         switch(it->first)
         {
            case BoundaryConditions::VALUE:
               lines.real().col(idx) = bcOp.value(it->second);
               idx++;
               hasReal = true;
               break;
            case BoundaryConditions::FIRST_DERIVATIVE:
               lines.real().col(idx) = bcOp.firstDerivative(it->second);
               idx++;
               hasReal = true;
               break;
            case BoundaryConditions::SECOND_DERIVATIVE:
               lines.real().col(idx) = bcOp.secondDerivative(it->second);
               idx++;
               hasReal = true;
               break;
            case BoundaryConditions::BETA_SLOPE:
               /// \warning Beta slope boundary conditions does not include the wave number factor
               if(it->second == IBoundary::RIGHT)
               {
                  lines.imag().col(idx) = static_cast<MHDFloat>(-1)*bcOp.value(it->second);
               } else //if(it->second == IBoundary::LEFT)
               {
                  lines.imag().col(idx) = bcOp.value(it->second);
               }
               idx++;
               hasImag = true;
               break;
            default:
               throw Exception("Unknown boundary condition ID for tau method");
               break;
         }
      }

      // Clear real matrix if not conditions have been set
      if(!hasReal)
      {
         lines.real().resize(0,0);
      }

      // Clear imaginary matrix if not conditions have been set
      if(!hasImag)
      {
         lines.imag().resize(0,0);
      }

      return lines;
   }

   DecoupledZSparse BoundaryConditions::tauMatrix(const IBoundary& bcOp, const std::vector<std::pair<BoundaryConditions::Id,IBoundary::Position> >& bcId)
   {
      // Builde the tau lines
      DecoupledZMatrix lines = BoundaryConditions::tauLines(bcOp, bcId);

      // Create sparse matrices
      DecoupledZSparse sparse(bcOp.basisN(), bcOp.basisN());

      // Reserve space in sparse matrices
      sparse.real().reserve(lines.real().rows()*lines.real().cols());
      sparse.imag().reserve(lines.real().rows()*lines.real().cols());

      // Make sparse tau matrix for real component (if required)
      if(sparse.real().size() != 0)
      {
         for(int j = 0; j < sparse.real().cols(); ++j)
         {
            // create column j
            sparse.real().startVec(j);

            for(int i = 0; i < lines.real().cols(); ++i)
            {
               sparse.real().insertBack(i,j) = lines.real()(j,i);
            }
         }
         sparse.real().prune(1e-32);
         sparse.real().finalize(); 
      }

      // Make sparse tau matrix for imaginary component (if required)
      if(sparse.imag().size() != 0)
      {
         for(int j = 0; j < sparse.imag().cols(); ++j)
         {
            // create column j
            sparse.imag().startVec(j);

            for(int i = 0; i < lines.imag().cols(); ++i)
            {
               sparse.imag().insertBack(i,j) = lines.imag()(j,i);
            }
         }
         sparse.imag().prune(1e-32);
         sparse.imag().finalize(); 
      }

      return sparse;
   }

}
}
