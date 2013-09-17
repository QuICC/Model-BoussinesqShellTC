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
      DecoupledZMatrix lines(std::make_pair(Matrix(bcOp.basisN(), nBCs),Matrix(bcOp.basisN(), nBCs)));
      lines.first.setZero();
      lines.second.setZero();

      // Map iterator
      std::vector<std::pair<BoundaryConditions::Id,IBoundary::Position> >::const_iterator it;

      // Create boundary values
      int idx = 0;
      for(it = bcId.begin(); it != bcId.end(); ++it)
      {
         switch(it->first)
         {
            case BoundaryConditions::VALUE:
               lines.first.col(idx) = bcOp.value(it->second);
               idx++;
               hasReal = true;
               break;
            case BoundaryConditions::FIRST_DERIVATIVE:
               lines.first.col(idx) = bcOp.firstDerivative(it->second);
               idx++;
               hasReal = true;
               break;
            case BoundaryConditions::SECOND_DERIVATIVE:
               lines.first.col(idx) = bcOp.secondDerivative(it->second);
               idx++;
               hasReal = true;
               break;
            case BoundaryConditions::BETA_SLOPE:
               /// \warning Beta slope boundary conditions does not include the wave number factor
               if(it->second == IBoundary::RIGHT)
               {
                  lines.second.col(idx) = static_cast<MHDFloat>(-1)*bcOp.value(it->second);
               } else //if(it->second == IBoundary::LEFT)
               {
                  lines.second.col(idx) = bcOp.value(it->second);
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
         lines.first.resize(0,0);
      }

      // Clear imaginary matrix if not conditions have been set
      if(!hasImag)
      {
         lines.second.resize(0,0);
      }

      return lines;
   }

   DecoupledZSparse BoundaryConditions::tauMatrix(const IBoundary& bcOp, const std::vector<std::pair<BoundaryConditions::Id,IBoundary::Position> >& bcId)
   {
      // Builde the tau lines
      DecoupledZMatrix lines = BoundaryConditions::tauLines(bcOp, bcId);

      // Create sparse matrices
      DecoupledZSparse sparse(std::make_pair(SparseMatrix(bcOp.basisN(), bcOp.basisN()),SparseMatrix(bcOp.basisN(), bcOp.basisN())));

      // Reserve space in sparse matrices
      sparse.first.reserve(lines.first.rows()*lines.first.cols());
      sparse.second.reserve(lines.first.rows()*lines.first.cols());

      // Make sparse tau matrix for real component (if required)
      if(sparse.first.size() != 0)
      {
         for(int j = 0; j < sparse.first.cols(); ++j)
         {
            // create column j
            sparse.first.startVec(j);

            for(int i = 0; i < lines.first.cols(); ++i)
            {
               sparse.first.insertBack(i,j) = lines.first(j,i);
            }
         }
         sparse.first.prune(1e-32);
         sparse.first.finalize(); 
      }

      // Make sparse tau matrix for imaginary component (if required)
      if(sparse.second.size() != 0)
      {
         for(int j = 0; j < sparse.second.cols(); ++j)
         {
            // create column j
            sparse.second.startVec(j);

            for(int i = 0; i < lines.second.cols(); ++i)
            {
               sparse.second.insertBack(i,j) = lines.second(j,i);
            }
         }
         sparse.second.prune(1e-32);
         sparse.second.finalize(); 
      }

      return sparse;
   }

}
}
