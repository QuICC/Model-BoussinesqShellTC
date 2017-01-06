/** 
 * @file SparseSolverDataSelector.hpp
 * @brief Small template to select the right data type for the sparse solvers
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSESOLVERDATASELECTOR_HPP
#define SPARSESOLVERDATASELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "SparseSolvers/SparseDummySolver.hpp"

namespace QuICC {

   namespace Solver {

      // Configure code to use TTT scheme
      #if defined QUICC_SPATIALSCHEME_TTT || defined QUICC_SPATIALSCHEME_TT
         template <template <class,class> class TSolver> class SparseSolverDataSelector
         {
            public:
               /// Typedef for a real operator solver
               typedef TSolver<SparseMatrix,Matrix>  RealSolverType;

               /// Typedef for a complex operator solver
               typedef SparseDummySolver<SparseMatrixZ,MatrixZ>  ComplexSolverType;
         };

      #else
         template <template <class,class> class TSolver> class SparseSolverDataSelector
         {
            public:
               /// Typedef for a real operator solver
               typedef TSolver<SparseMatrix,DecoupledZMatrix>  RealSolverType;

               /// Typedef for a complex operator solver
               typedef TSolver<SparseMatrixZ,MatrixZ>  ComplexSolverType;
         };
      #endif //QUICC_SPATIALSCHEME_TTT

   }
}

#endif // SPARSESOLVERDATASELECTOR_HPP
