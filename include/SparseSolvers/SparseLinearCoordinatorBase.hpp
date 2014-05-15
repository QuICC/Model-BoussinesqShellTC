/** 
 * @file SparseLinearCoordinatorBase.hpp
 * @brief Implementation of the base for a general sparse linear solver coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSELINEARCOORDINATORBASE_HPP
#define SPARSELINEARCOORDINATORBASE_HPP

// Debug includes
//
#include "Debug/DebuggerMacro.h"

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Project includes
//
#include "Base/MathConstants.hpp"
#include "SparseSolvers/SparseCoordinatorBase.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Solver {

   namespace internal {

      void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat);

      void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat);
   }

   /**
    * @brief Implementation of the base for a general sparse linear solver coordinator
    */
   template <template <class,class> class TSolver, bool TIsComplex> class SparseLinearCoordinatorBase;

   namespace internal {

      inline void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat)
      {
         assert(decMat.real().size() > 0);
         assert(decMat.imag().size() == 0 || decMat.imag().nonZeros() == 0);

         if(c != 1.0)
         {
            mat += c*decMat.real();
         } else
         {
            mat += decMat.real();
         }
      }

      inline void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat)
      {
         assert(decMat.real().size() > 0);
         assert(decMat.imag().size() > 0);
         assert(decMat.real().size() == decMat.imag().size());

         if(c != 1.0)
         {
            mat += c*decMat.real().cast<MHDComplex>() + c*Math::cI*decMat.imag();
         } else
         {
            mat += decMat.real().cast<MHDComplex>() + Math::cI*decMat.imag();
         }
      }
   }
}
}

// Include complex field spezialization
#include "SparseSolvers/SparseLinearZCoordinatorBase.hpp"

// Include real field spezialization
#include "SparseSolvers/SparseLinearRCoordinatorBase.hpp"

#endif // SPARSELINEARCOORDINATORBASE_HPP
