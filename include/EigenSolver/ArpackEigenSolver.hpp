/** \file ArpackEigenSolver.hpp
 *  \brief Definition of an interface to the ARPACK sparse eigen solver
 */

#ifndef ARPACKEIGENSOLVER_HPP
#define ARPACKEIGENSOLVER_HPP

// System includes
//
#include <cmath>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Timesteppers/SparseSolverMacro.h"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of an eigen solver based on ARPACK
    */
   class ArpackEigenSolver
   {
      public:
         /**
          * @brief constructor
          */
         ArpackEigenSolver();

         /**
          * @brief Destructor
          */
         ~ArpackEigenSolver();

         /**
          * @brief Setup and factorise the LHS and RHS matrices of the eigenvalue problem
          */
         void compute(const SparseMatrixZ& lhs, const SparseMatrixZ& rhs);

         /**
          * @brief Solve the eigen problem
          */
         void solve(ArrayZ& eigenValues, MatrixZ& eigenVectors);

         /**
          * @brief Get the status from ARPACK
          */
         int info() const;

         /**
          * @brief Set the shift eigenvalue
          */
         void setSigma(const MHDComplex sigma);

         /**
          * @brief Set the tolerance
          */
         void setTol(const MHDFloat tol);

         /**
          * @brief Set the maximum number of Arnoldi iterations
          */
         void setMaxIter(const int maxIter);

         /**
          * @brief Set NCV
          */
         void setNcv(const int ncv);

         /**
          * @brief Set the type of eigenvalues to look for
          */
         void setWhich(const std::string& which);
         
      protected:

      private:
         /**
          * @brief Allocate work memory
          */
         void allocateMemory();

         /**
          * @brief Storage for the status output from ARPACK
          */
         int mInfo;

         /**
          * @brief Shift eigenvalue \f$\sigma\f$
          */
         MHDComplex  mSigma;

         /**
          * @brief TOL option for ARPACK
          */
         MHDFloat mTol;

         /**
          * @brief MXITER option for ARPACK
          */
         int mMaxIter;

         /**
          * @brief NCV option for ARPACK
          */
         int mNcv;

         /**
          * @brief WHICH option for ARPACK
          */
         std::string mWhich;

         /**
          * @brief LHS operator of the eigenvalue problem
          */
         SparseMatrixZ mLhs;

         /**
          * @brief RHS operator of the eigenvalue problem
          */
         SparseMatrixZ mRhs;

         /**
          * @brief Sparse linear solver
          */
         SparseSolverMacro<SparseMatrixZ> mSolver;

         /**
          * @defgroup Work memory for ARPACK
          */
         /**@{*/
         ArrayI mIpntr;
         ArrayI mIparam;
         ArrayZ mResid;
         MatrixZ mV;
         ArrayZ mWorkd;
         ArrayZ mWorkl;
         Array mRwork;
         ArrayI mSelect;
         ArrayZ mD;
         MatrixZ mZ;
         ArrayZ mWorkev;
         /**@}*/

   };

   inline bool sortEigenValues(const MHDComplex& x, const MHDComplex& y)
   {
      if(std::isinf(x.real()) || std::abs(x.real()) > 1e20)
      {
         return false;
      }
      else if(std::isinf(y.real()) || std::abs(x.real()) > 1e20)
      {
         return true;
      } else
      {
         return x.real() > y.real();
      }
   }
}

#endif // ARPACKEIGENSOLVER_HPP
