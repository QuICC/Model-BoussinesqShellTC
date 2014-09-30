/**
 * @file MumpsLU.hpp
 * @brief Eigen wrapper to use the Mumps solver 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MUMPSLU_HPP
#define MUMPSLU_HPP

#include "Eigen/SparseCore"
#include "Eigen/src/Core/util/DisableStupidWarnings.h"

#include "Eigen/src/misc/Solve.h"
#include "Eigen/src/misc/SparseSolve.h"

#include <complex>
#include <smumps_c.h>
#include <cmumps_c.h>
#include <dmumps_c.h>
#include <zmumps_c.h>

#include "Framework/FrameworkMacro.h"
#include "Base/MpiTypes.hpp"

namespace Eigen {

   #define MUMPS_FORTRAN_COMM -987654

   template <typename TMatrixType> class MumpsStrucType;
   template <typename TScalar> class MumpsScalarType
   {
      public:
         typedef TScalar ScalarType;
   };
   template <> class MumpsScalarType<std::complex<float> >
   {
      public:
         typedef mumps_complex ScalarType;
   };
   template <> class MumpsScalarType<std::complex<double> >
   {
      public:
         typedef mumps_double_complex ScalarType;
   };

   #define DECL_MUMPS(PREFIX,UPREFIX,KEYTYPE)		                                             \
      template <> class MumpsStrucType<KEYTYPE>                                                 \
      {                                                                                         \
         public:                                                                                \
            typedef UPREFIX##MUMPS_STRUC_C StrucType;                                           \
      };                                                                                        \
      extern "C" {                                                                              \
         extern void PREFIX##mumps_c(MumpsStrucType<KEYTYPE>::StrucType * idptr);               \
      }                                                                                         \
      inline void MumpsLU_mumps(MumpsStrucType<KEYTYPE>::StrucType * idptr, KEYTYPE)            \
      {                                                                                         \
         PREFIX##mumps_c(idptr);                                                                \
      }                                                                                         \

   DECL_MUMPS(s,S,float)
   DECL_MUMPS(c,C,std::complex<float>)
   DECL_MUMPS(d,D,double)
   DECL_MUMPS(z,Z,std::complex<double>)

   inline float MumpsLU_convert(const float val) {return val;};
   inline double MumpsLU_convert(const double val) {return val;};
   inline mumps_complex MumpsLU_convert(const std::complex<float> val) {return (mumps_complex){.r = val.real(), .i = val.imag()};};
   inline mumps_double_complex MumpsLU_convert(const std::complex<double> val) {return (mumps_double_complex){.r = val.real(), .i = val.imag()};};
   inline std::complex<float> MumpsLU_convert(const mumps_complex val) {return std::complex<float>(val.r, val.i);};
   inline std::complex<double> MumpsLU_convert(const mumps_double_complex val) {return std::complex<float>(val.r, val.i);};

   template <typename TMatrixType>
   class MumpsLU: internal::noncopyable
   {
      public:
         typedef TMatrixType MatrixType;
         typedef typename MatrixType::Scalar Scalar;
         typedef typename MatrixType::RealScalar RealScalar;
         typedef typename MatrixType::Index Index;
         typedef Matrix<Scalar,Dynamic,1> Vector;
         typedef Matrix<int, 1, MatrixType::ColsAtCompileTime> IntRowVectorType;
         typedef Matrix<int, MatrixType::RowsAtCompileTime, 1> IntColVectorType;
         typedef typename MumpsScalarType<Scalar>::ScalarType MumpsScalar;
         typedef SparseMatrix<Scalar,RowMajor,int> MumpsMatrixType;

      public:

         MumpsLU() { init(); }

         MumpsLU(const MatrixType& matrix)
         {
            init();
            compute(matrix);
         }

         ~MumpsLU()
         {
            GeoMHDiSCC::FrameworkMacro::synchronize();
            m_id.job = -2;

            #ifdef GEOMHDISCC_MPI
               free(m_id.irn_loc);
               free(m_id.jcn_loc);
               free(m_id.a_loc);
            #else 
               free(m_id.irn);
               free(m_id.jcn);
               free(m_id.a);
            #endif //GEOMHDISCC_MPI
            free(m_id.rhs);

            MumpsLU_mumps(&m_id, Scalar());
            GeoMHDiSCC::FrameworkMacro::synchronize();
         }

         inline Index rows() const { return m_id.n; }
         inline Index cols() const { return m_id.n; }

         inline void setIterativeRefinement(const int iterations) { m_id.icntl[10 - 1] = iterations; }

         inline void setMemoryRelaxation(const int percent) { m_id.icntl[14 - 1] = percent; }

         inline int getMaxMemory() { return m_id.infog[16 - 1]; }
         inline int getTotalMemory() { return m_id.infog[17 - 1]; }

         /** \brief Reports whether previous computation was successful.
          *
          * \returns \c Success if computation was succesful,
          *          \c NumericalIssue if the matrix.appears to be negative.
          */
         ComputationInfo info() const
         {
            eigen_assert(m_isInitialized && "Decomposition is not initialized.");
            return m_info;
         }

         /** Computes the sparse LU decomposition of \a matrix 
          *  Note that the matrix should be column-major, and in compressed format for best performance.
          *  \sa SparseMatrix::makeCompressed().
          */
         void compute(const MatrixType& matrix)
         {
            analyzePattern(matrix);
            factorize(matrix);
         }

         /** \returns the solution x of \f$ A x = b \f$ using the current decomposition of A.
          *
          * \sa compute()
          */
         template<typename Rhs> inline const internal::solve_retval<MumpsLU, Rhs> solve(const MatrixBase<Rhs>& b) const
         {
            eigen_assert(m_isInitialized && "MumpsLU is not initialized.");
            eigen_assert(rows()==b.rows()
                  && "MumpsLU::solve(): invalid number of rows of the right hand side matrix b");
            return internal::solve_retval<MumpsLU, Rhs>(*this, b.derived());
         }

         /** \returns the solution x of \f$ A x = b \f$ using the current decomposition of A.
          *
          * \sa compute()
          */
         template<typename Rhs> inline const internal::sparse_solve_retval<MumpsLU, Rhs> solve(const SparseMatrixBase<Rhs>& b) const
         {
            eigen_assert(m_isInitialized && "MumpsLU is not initialized.");
            eigen_assert(rows()==b.rows()
                  && "MumpsLU::solve(): invalid number of rows of the right hand side matrix b");
            return internal::sparse_solve_retval<MumpsLU, Rhs>(*this, b.derived());
         }

         /** Performs a symbolic decomposition on the sparcity of \a matrix.
          *
          * This function is particularly useful when solving for several problems having the same structure.
          *
          * \sa factorize(), compute()
          */
         void analyzePattern(const MatrixType& matrix)
         {
            GeoMHDiSCC::FrameworkMacro::synchronize();

            m_id.job = 1;
            copyInput(matrix, true);

            MumpsLU_mumps(&m_id, Scalar());
            m_error = m_id.info[1-1];

            m_isInitialized = true;
            m_info = m_error ? InvalidInput : Success;
            m_analysisIsOk = m_error ? false : true;
            m_factorizationIsOk = false;
            GeoMHDiSCC::FrameworkMacro::synchronize();
         }

         /** Performs a numeric decomposition of \a matrix
          *
          * The given matrix must have the same sparcity than the matrix on which the pattern analysis has been performed.
          *
          * \sa analyzePattern(), compute()
          */
         void factorize(const MatrixType& matrix)
         {
            GeoMHDiSCC::FrameworkMacro::synchronize();

            eigen_assert(m_analysisIsOk && "MumpsLU: you must first call analyzePattern()");

            m_id.job = 2;
            copyInput(matrix, false);

            MumpsLU_mumps(&m_id, Scalar());
            m_error = m_id.info[1-1];

            m_info = m_error ? NumericalIssue : Success;
            m_factorizationIsOk = m_error ? false : true;
            GeoMHDiSCC::FrameworkMacro::synchronize();
         }

         #ifndef EIGEN_PARSED_BY_DOXYGEN
         /** \internal */
         template<typename BDerived,typename XDerived> bool _solve(const MatrixBase<BDerived> &b, MatrixBase<XDerived> &x) const;
         #endif

      protected:

         void init()
         {
            GeoMHDiSCC::FrameworkMacro::synchronize();

            m_id.job = -1;
            #ifdef GEOMHDISCC_MPI
               m_id.comm_fortran = MPI_Comm_c2f(GeoMHDiSCC::FrameworkMacro::spectralComm());
               int rank;
               MPI_Comm_rank(GeoMHDiSCC::FrameworkMacro::spectralComm(), &rank);
               m_isHost = (rank == 0);
            #else
               m_id.comm_fortran = MUMPS_FORTRAN_COMM;
               m_isHost = 1;
            #endif // GEOMHDISCC_MPI
            m_id.par = 1;
            m_id.sym = 0;

            m_info = InvalidInput;
            m_isInitialized  = false;
            m_error          = 0;

            MumpsLU_mumps(&m_id, Scalar());
            m_error = m_id.info[1-1];
            #ifdef GEOMHDISCC_NO_DEBUG
               m_id.icntl[1-1] = 0;
               m_id.icntl[2-1] = 0;
               m_id.icntl[3-1] = 0;
               m_id.icntl[4-1] = 0;
            #endif // GEOMHDISCC_NO_DEBUG

            #ifdef GEOMHDISCC_MPI
               // Matrix distribution
               m_id.icntl[18-1] = 3;
               // RHS dense
               m_id.icntl[20-1] = 0;
               // Solution distribution
               m_id.icntl[21-1] = 0;
            #endif // GEOMHDISCC_MPI

            // increase memory relaxation
            m_id.icntl[14-1] = 35;

            GeoMHDiSCC::FrameworkMacro::synchronize();
         }

         void copyInput(const MatrixType& mat, const bool freeMemory)
         {
            eigen_assert(mat.rows() == mat.cols() && "MumpsLU requires a square matrix");

            if(m_id.a != 0 && freeMemory)
            {
               #ifdef GEOMHDISCC_MPI
                  free(m_id.irn_loc);
                  free(m_id.jcn_loc);
                  free(m_id.a_loc);
                  m_id.irn_loc = 0;
                  m_id.jcn_loc = 0;
                  m_id.a_loc = 0;
               #else
                  free(m_id.irn);
                  free(m_id.jcn);
                  free(m_id.a);
                  m_id.irn = 0;
                  m_id.jcn = 0;
                  m_id.a = 0;
               #endif //GEOMHDISCC_MPI

               free(m_id.rhs);
               m_id.rhs = 0;
            }

            bool needCopy;
            #ifdef GEOMHDISCC_MPI
               needCopy = (m_id.a_loc == 0);
            #else 
               needCopy = (m_id.a == 0);
            #endif //GEOMHDISCC_MPI

            if(needCopy)
            {
               int *pRow;
               int *pCol;
               MumpsScalar *pValue;

               m_id.n = mat.rows();
               #ifdef GEOMHDISCC_MPI
                  m_id.nz_loc = mat.nonZeros();
                  m_id.irn_loc = (int *) malloc(sizeof(int)*m_id.nz_loc);
                  m_id.jcn_loc = (int *) malloc(sizeof(int)*m_id.nz_loc);

                  // Mumps requires a different format for complex values
                  m_id.a_loc = (MumpsScalar *) malloc(sizeof(MumpsScalar)*m_id.nz_loc);

                  // Set pointers to copy data
                  pRow = m_id.irn_loc;
                  pCol = m_id.jcn_loc;
                  pValue = m_id.a_loc;
               #else
                  m_id.nz = mat.nonZeros();
                  m_id.irn = (int *) malloc(sizeof(int)*m_id.nz);
                  m_id.jcn = (int *) malloc(sizeof(int)*m_id.nz);

                  // Mumps requires a different format for complex values
                  m_id.a = (MumpsScalar *) malloc(sizeof(MumpsScalar)*m_id.nz);

                  // Set pointers to copy data
                  pRow = m_id.irn;
                  pCol = m_id.jcn;
                  pValue = m_id.a;
               #endif //GEOMHDISCC_MPI

               // Also allocate memory for rhs on host, assumes single rhs
               if(m_isHost)
               {
                  m_id.rhs = (MumpsScalar *) malloc(sizeof(MumpsScalar)*m_id.n);
               }

               for(int k = 0; k < mat.outerSize(); ++k)
               {
                  for(typename MatrixType::InnerIterator it(mat,k); it; ++it)
                  {
                     // Convert index to fortran numbering (1 based vs 0 based)
                     *pRow = it.row() + 1;
                     *pCol = it.col() + 1;
                     // Convert value to compatible data type
                     *pValue = MumpsLU_convert(it.value());

                     // Increment pointers
                     ++pRow;
                     ++pCol;
                     ++pValue;
                  }
               }
            }
         }

         void copyRhs(const Scalar* pData) const
         {
            const Scalar* pRhs;

            #ifdef GEOMHDISCC_MPI
               mTmp.resize(m_id.lrhs, m_id.nrhs);
               MPI_Reduce(pData, mTmp.data(), m_id.nrhs*m_id.lrhs, GeoMHDiSCC::Parallel::MpiTypes::type<Scalar>(), MPI_SUM, 0, GeoMHDiSCC::FrameworkMacro::spectralComm()); 
               pRhs = mTmp.data();
            #else
               pRhs = pData;
            #endif //GEOMHDISCC_MPI

            if(m_isHost)
            {
               int nK = m_id.nrhs*m_id.lrhs;
               MumpsScalar * pVal = m_id.rhs;

               for(int k = 0; k < nK; ++pRhs,++pVal,++k)
               {
                  *pVal = MumpsLU_convert(*pRhs);
               }
            }
         }

         void copySolution(Scalar* pData) const
         {
            int nK = m_id.nrhs*m_id.lrhs;
            if(m_isHost)
            {
               Scalar* pSol = pData;
               MumpsScalar * pVal = m_id.rhs;

               for(int k = 0; k < nK; ++pSol,++pVal,++k)
               {
                  *pSol = MumpsLU_convert(*pVal);
               }
            }

            #ifdef GEOMHDISCC_MPI
               MPI_Bcast(pData, nK, GeoMHDiSCC::Parallel::MpiTypes::type<Scalar>(), 0, GeoMHDiSCC::FrameworkMacro::spectralComm());
            #endif //GEOMHDISCC_MPI
         }

         mutable ComputationInfo m_info;
         bool m_isInitialized;
         int m_factorizationIsOk;
         int m_analysisIsOk;

         mutable int m_error;
         mutable typename MumpsStrucType<Scalar>::StrucType m_id;

         int m_isHost;

      private:
         mutable Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>   mTmp;

         MumpsLU(MumpsLU& ) { }

   };

template<typename MatrixType>
template<typename BDerived,typename XDerived>
bool MumpsLU<MatrixType>::_solve(const MatrixBase<BDerived> &b, MatrixBase<XDerived> &x) const
{
   GeoMHDiSCC::FrameworkMacro::synchronize();

   int rhsCols = b.cols();
   eigen_assert((BDerived::Flags&RowMajorBit)==0 && "MumpsLU backend does not support non col-major rhs yet");
   eigen_assert((XDerived::Flags&RowMajorBit)==0 && "MumpsLU backend does not support non col-major result yet");
   eigen_assert(b.cols() == 1 && "MumpsLU is implemented for a single RHS");
   #ifdef GEOMHDISCC_MPI
      eigen_assert(b.derived().data() != x.derived().data() && "Parallel execution does not allow x and b to be the same");
   #endif //GEOMHDISCC_MPI

   m_id.job = 3;
   m_id.nrhs = rhsCols;
   m_id.lrhs = m_id.n;
   m_id.icntl[10 - 1] = 3;
   #ifdef GEOMHDISCC_DEBUG
      m_id.icntl[11-1] = 1;
   #endif // GEOMHDISCC_DEBUG

   copyRhs(b.derived().data());

   MumpsLU_mumps(&m_id, Scalar());
   m_error = m_id.info[1-1];

   if (m_error!=0)
      return false;

   copySolution(x.derived().data());

   GeoMHDiSCC::FrameworkMacro::synchronize();

   return true;
}

namespace internal {

template<typename TMatrixType, typename Rhs>
struct solve_retval<MumpsLU<TMatrixType>, Rhs>
  : solve_retval_base<MumpsLU<TMatrixType>, Rhs>
{
  typedef MumpsLU<TMatrixType> Dec;
  EIGEN_MAKE_SOLVE_HELPERS(Dec,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    dec()._solve(rhs(),dst);
  }
};

template<typename TMatrixType, typename Rhs>
struct sparse_solve_retval<MumpsLU<TMatrixType>, Rhs>
  : sparse_solve_retval_base<MumpsLU<TMatrixType>, Rhs>
{
  typedef MumpsLU<TMatrixType> Dec;
  EIGEN_MAKE_SPARSE_SOLVE_HELPERS(Dec,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    this->defaultEvalTo(dst);
  }
};

} // end namespace internal

} // end namespace Eigen

#include "Eigen/src/Core/util/ReenableStupidWarnings.h"

#endif // MUMPSLU_HPP
