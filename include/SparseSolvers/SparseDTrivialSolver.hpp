/** 
 * @file SparseDTrivialSolver.hpp
 * @brief Implementation of a real (coupled) trivial solver structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef SPARSEDTRIVIALSOLVER_HPP
#define SPARSEDTRIVIALSOLVER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "SparseSolvers/SparseSolverBase.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of a real (coupled) trivial solver structure
    */
   class SparseDTrivialSolver: public SparseSolverBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          */
         SparseDTrivialSolver(const int start);

         /**
          * @brief Destructor
          */
         virtual ~SparseDTrivialSolver();

         /**
          * @brief Get the number of systems in solver
          */
         int nSystem() const;

         /**
          * @brief Add data storage
          * 
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void addStorage(const int rows, const int cols);

         /**
          * @brief Set RHS data
          *
          * @param idx   Index of the data
          */
         DecoupledZMatrix& rRHSData(const int idx);

         /**
          * @brief Get solution data
          *
          * @param idx   Index of the data
          */
         const DecoupledZMatrix& solution(const int idx) const;

         /**
          * @brief Set solution data
          *
          * @param idx   Index of the data
          */
         DecoupledZMatrix& rSolution(const int idx);
         
      protected:
         /**
          * @brief Storage for solution of linear solve
          */
         std::vector<DecoupledZMatrix>  mSolution;

      private:
   };

   /// Typedef for a shared pointer of a SparseDTrivialSolver
   typedef SharedPtrMacro<SparseDTrivialSolver>  SharedSparseDTrivialSolver;
}
}

#endif // SPARSEDTRIVIALSOLVER_HPP
