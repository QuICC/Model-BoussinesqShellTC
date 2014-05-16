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

namespace GeoMHDiSCC {

   namespace Solver {

      // Configure code to use TTT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TTT
         template <template <class,class> class TSolver> class SparseSolverDataSelector
         {
            public:
               /// Typedef for a real operator solver
               typedef TSolver<SparseMatrix,Matrix>  RealSolverType;

               /// Typedef for a complex operator solver
               typedef DummySolver<SparseMatrixZ,MatrixZ>  ComplexSolverType;
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
      #endif //GEOMHDISCC_SPATIALSCHEME_TTT

   }
}

#endif // SPARSESOLVERDATASELECTOR_HPP
