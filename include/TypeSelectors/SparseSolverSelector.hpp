/** \file SparseSolverSelector.hpp
 *  \brief Small template to select different sparse solvers.
 */

#ifndef SPARSESOLVERSELECTOR_HPP
#define SPARSESOLVERSELECTOR_HPP

// SuperLU sparse solver
#ifdef GEOMHDISCC_SPLINALG_SUPERLU
   // Include the right header
   #include <Eigen/SuperLUSupport>

   namespace GeoMHDiSCC {
      /**
       * @brief Selector for the different implementations of the sparse solvers. Sets up the SuperLU version. 
       */
      template <typename TMatrix> struct SparseSolverSelector
      {
         typedef Eigen::SuperLU<TMatrix> SolverType;
      };
   }
#endif //GEOMHDISCC_SPLINALG_SUPERLU

// UmfPack Version
#ifdef GEOMHDISCC_SPLINALG_UMFPACK
   // Include the right header
   #include <Eigen/UmfPackSupport>

   namespace GeoMHDiSCC {
      /**
       * @brief Selector for the different implementations of the sparse solvers. Sets up the UmfPack version. 
       */
      template <typename TMatrix> struct SparseSolverSelector
      {
         typedef Eigen::UmfPackLU<TMatrix> SolverType;
      };
   }
#endif //GEOMHDISCC_SPLINALG_UMFPACK

// SparseLU Version
#ifdef GEOMHDISCC_SPLINALG_SPARSELU
   // Include the right header
   #include <Eigen/SparseLU>
   #include <Eigen/OrderingMethods>
   //#include <Eigen/MetisSupport>

   namespace GeoMHDiSCC {
      /**
       * @brief Selector for the different implementations of the sparse solvers. Sets up the UmfPack version. 
       */
      template <typename TMatrix> struct SparseSolverSelector
      {
         //typedef Eigen::SparseLU<TMatrix, Eigen::NaturalOrdering<typename TMatrix::Index> > SolverType;
         //typedef Eigen::SparseLU<TMatrix, Eigen::AMDOrdering<typename TMatrix::Index> > SolverType;
         typedef Eigen::SparseLU<TMatrix, Eigen::COLAMDOrdering<typename TMatrix::Index> > SolverType;
         //typedef Eigen::SparseLU<TMatrix, Eigen::MetisOrdering<typename TMatrix::Index> > SolverType;
      };
   }
#endif //GEOMHDISCC_SPLINALG_SPARSELU

// KLU Version
#ifdef GEOMHDISCC_SPLINALG_KENTLU
   // Include the right header
   #include <Eigen/KentLUSupport>

   namespace GeoMHDiSCC {
      /**
       * @brief Selector for the different implementations of the sparse solvers. Sets up the KentLU version. 
       */
      template <typename TMatrix> struct SparseSolverSelector
      {
         typedef Eigen::KentLU<TMatrix> SolverType;
      };
   }
#endif //GEOMHDISCC_SPLINALG_KENTLU

// SparseQR Version
#ifdef GEOMHDISCC_SPLINALG_SPARSEQR
   // Include the right header
   #include <Eigen/SparseQR>
   #include <Eigen/OrderingMethods>
   //#include <Eigen/MetisSupport>

   namespace GeoMHDiSCC {
      /**
       * @brief Selector for the different implementations of the sparse solvers. Sets up the SparseQR version. 
       */
      template<typename TMatrix> struct SparseSolverSelector
      {
         //typedef Eigen::SparseQR<TMatrix, Eigen::NaturalOrdering<typename TMatrix::Index> > SolverType;
         //typedef Eigen::SparseQR<TMatrix, Eigen::AMDOrdering<typename TMatrix::Index> > SolverType;
         typedef Eigen::SparseQR<TMatrix, Eigen::COLAMDOrdering<typename TMatrix::Index> > SolverType;
         //typedef Eigen::SparseQR<TMatrix, Eigen::MetisOrdering<typename TMatrix::Index> > SolverType;
      };
   }
#endif //GEOMHDISCC_SPLINALG_SPARSEQR

// SuiteSparseQR Version
#ifdef GEOMHDISCC_SPLINALG_SPQR
   // Include the right header
   #include <Eigen/SPQRSupport>

   namespace GeoMHDiSCC {
      /**
       * @brief Selector for the different implementations of the sparse solvers. Sets up the SuiteSparseQR version. 
       */
      template<typename TMatrix> struct SparseSolverSelector
      {
         typedef Eigen::SPQR<TMatrix> SolverType;
      };
   }
#endif //GEOMHDISCC_SPLINALG_SPQR

// Pardiso Version
#ifdef GEOMHDISCC_SPLINALG_PARDISO
   // Include the right header
   #include "../External/Interfaces/PardisoLU.hpp"

   namespace GeoMHDiSCC {
      /**
       * @brief Selector for the different implementations of the sparse solvers. Sets up the Pardiso version. 
       */
      template<typename TMatrix> struct SparseSolverSelector
      {
         typedef Eigen::PardisoLU<TMatrix> SolverType;
      };
   }
#endif //GEOMHDISCC_SPLINALG_PARDISO

// MKL Pardiso Version
#ifdef GEOMHDISCC_SPLINALG_MKLPARDISO
   // Include the right header
   #include <Eigen/PardisoSupport>

   namespace GeoMHDiSCC {
      /**
       * @brief Selector for the different implementations of the sparse solvers. Sets up the Pardiso version. 
       */
      template<typename TMatrix> struct SparseSolverSelector
      {
         typedef Eigen::PardisoLU<TMatrix> SolverType;
      };
   }
#endif //GEOMHDISCC_SPLINALG_MKLPARDISO

// BiCGSTAB Version
#ifdef GEOMHDISCC_SPLINALG_BICGSTAB
   // Include the right header
   #include <Eigen/IterativeLinearSolvers>   

   namespace GeoMHDiSCC {
      /**
       * @brief Selector for the different implementations of the sparse solvers. Sets up the BiCGSTAB version. 
       */
      template<typename TMatrix> struct SparseSolverSelector
      {
         //typedef Eigen::BiCGSTAB<TMatrix, Eigen::IdentityPreconditioner<typename TMatrix::Scalar> > SolverType;
         typedef Eigen::BiCGSTAB<TMatrix, Eigen::DiagonalPreconditioner<typename TMatrix::Scalar> > SolverType;
         //typedef Eigen::BiCGSTAB<TMatrix, Eigen::IncompleteLUT<typename TMatrix::Scalar> > SolverType;
      };
   }
#endif //GEOMHDISCC_SPLINALG_BICGSTAB

#endif // SPARSESOLVERSELECTOR_HPP
