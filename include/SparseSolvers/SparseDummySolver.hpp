/** 
 * @file SparseDummySolver.hpp
 * @brief Implementation of the sparse dummy solver
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSEDUMMYSOLVER_HPP
#define SPARSEDUMMYSOLVER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/SolveTiming.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   /**
    * @brief Implementation of the sparse dummy solver
    */
   template <typename TOperator, typename TData> class SparseDummySolver
   {
      public:
         /**
          * @brief Dummy constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param time    Solver timing with respect to timestepping
          */
         SparseDummySolver(const int start, const SolveTiming::Id time) {};

         /**
          * @brief Dummy destructor
          */
         virtual ~SparseDummySolver() {};

         /**
          * @brief Dummy implementation 
          */
         void initStartRow() {};
         
      protected:

      private:
   };

   typedef SparseDummySolver<SparseMatrixZ,MatrixZ> SparseDummySolverComplexType;

   typedef std::vector<SharedPtrMacro<SparseDummySolverComplexType> >::iterator ComplexDummy_iterator;
}
}

#endif // SPARSEDUMMYSOLVER_HPP
