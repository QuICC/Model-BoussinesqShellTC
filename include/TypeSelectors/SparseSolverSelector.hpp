/** 
 * @file SparseSolverSelector.hpp
 * @brief Small template to select different sparse solvers.
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPARSESOLVERSELECTOR_HPP
#define SPARSESOLVERSELECTOR_HPP

// SuperLU Version for general solve
#ifdef GEOMHDISCC_SPLINALG_SUPERLU
   // Include the right header
   #include <Eigen/SuperLUSupport>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up the SuperLU version. 
          */
         template <typename TMatrix> struct SparseSelector
         {
            typedef Eigen::SuperLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_SUPERLU

// UmfPack Version for general solve
#ifdef GEOMHDISCC_SPLINALG_UMFPACK
   // Include the right header
   #include <Eigen/UmfPackSupport>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up the UmfPack version. 
          */
         template <typename TMatrix> struct SparseSelector
         {
            typedef Eigen::UmfPackLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_UMFPACK

// SparseLU Version for general solve
#ifdef GEOMHDISCC_SPLINALG_SPARSELU
   // Include the right header
   #include <Eigen/SparseLU>
   #include <Eigen/OrderingMethods>
   //#include <Eigen/MetisSupport>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up Eigen's SparseLU version. 
          */
         template <typename TMatrix> struct SparseSelector
         {
            //typedef Eigen::SparseLU<TMatrix, Eigen::NaturalOrdering<typename TMatrix::Index> > Type;
            //typedef Eigen::SparseLU<TMatrix, Eigen::AMDOrdering<typename TMatrix::Index> > Type;
            typedef Eigen::SparseLU<TMatrix, Eigen::COLAMDOrdering<typename TMatrix::Index> > Type;
            //typedef Eigen::SparseLU<TMatrix, Eigen::MetisOrdering<typename TMatrix::Index> > Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_SPARSELU

// KLU Version for general solve
#ifdef GEOMHDISCC_SPLINALG_KENTLU
   // Include the right header
   #include <Eigen/KentLUSupport>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up the KentLU version. 
          */
         template <typename TMatrix> struct SparseSelector
         {
            typedef Eigen::KentLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_KENTLU

// SparseQR Version for general solve
#ifdef GEOMHDISCC_SPLINALG_SPARSEQR
   // Include the right header
   #include <Eigen/SparseQR>
   #include <Eigen/OrderingMethods>
   //#include <Eigen/MetisSupport>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up the SparseQR version. 
          */
         template<typename TMatrix> struct SparseSelector
         {
            //typedef Eigen::SparseQR<TMatrix, Eigen::NaturalOrdering<typename TMatrix::Index> > Type;
            //typedef Eigen::SparseQR<TMatrix, Eigen::AMDOrdering<typename TMatrix::Index> > Type;
            typedef Eigen::SparseQR<TMatrix, Eigen::COLAMDOrdering<typename TMatrix::Index> > Type;
            //typedef Eigen::SparseQR<TMatrix, Eigen::MetisOrdering<typename TMatrix::Index> > Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_SPARSEQR

// SuiteSparseQR Version for general solve
#ifdef GEOMHDISCC_SPLINALG_SPQR
   // Include the right header
   #include <Eigen/SPQRSupport>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up the SuiteSparseQR version. 
          */
         template<typename TMatrix> struct SparseSelector
         {
            typedef Eigen::SPQR<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_SPQR

// Pardiso Version for general solve
#ifdef GEOMHDISCC_SPLINALG_PARDISO
   // Include the right header
   #include "../External/Interfaces/PardisoLU.hpp"

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up the Pardiso version. 
          */
         template<typename TMatrix> struct SparseSelector
         {
            typedef Eigen::PardisoLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_PARDISO

// MKL Pardiso Version for general solve
#ifdef GEOMHDISCC_SPLINALG_MKLPARDISO
   // Include the right header
   #include <Eigen/PardisoSupport>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up the MKL Pardiso version. 
          */
         template<typename TMatrix> struct SparseSelector
         {
            typedef Eigen::PardisoLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_MKLPARDISO

// MUMPS Version for general solve
#ifdef GEOMHDISCC_SPLINALG_MUMPS
   // Include the right header
   #include "../External/Interfaces/MumpsLU.hpp"

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up the MUMPS version. 
          */
         template<typename TMatrix> struct SparseSelector
         {
            typedef Eigen::MumpsLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_MUMPS

// BiCGSTAB Version for general solve
#ifdef GEOMHDISCC_SPLINALG_BICGSTAB
   // Include the right header
   #include <Eigen/IterativeLinearSolvers>   

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers. Sets up the BiCGSTAB version. 
          */
         template<typename TMatrix> struct SparseSelector
         {
            //typedef Eigen::BiCGSTAB<TMatrix, Eigen::IdentityPreconditioner<typename TMatrix::Scalar> > Type;
            typedef Eigen::BiCGSTAB<TMatrix, Eigen::DiagonalPreconditioner<typename TMatrix::Scalar> > Type;
            //typedef Eigen::BiCGSTAB<TMatrix, Eigen::IncompleteLUT<typename TMatrix::Scalar> > Type;
         };
      }
   }
#endif //GEOMHDISCC_SPLINALG_BICGSTAB

// UmfPack Version for SPD solve
#ifdef GEOMHDISCC_SPSPDLINALG_UMFPACK
   // Include the right header
   #include <Eigen/UmfPackSupport>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers for SPD matries. Sets up the UmfPack version. 
          */
         template <typename TMatrix> struct SparseSpdSelector
         {
            typedef Eigen::UmfPackLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPSPDLINALG_UMFPACK

// SparseLU Version for SPD solve
#ifdef GEOMHDISCC_SPSPDLINALG_SPARSELU
   // Include the right header
   #include <Eigen/SparseLU>
   #include <Eigen/OrderingMethods>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers for SPD matries. Sets up Eigen's SparseLU version. 
          */
         template <typename TMatrix> struct SparseSpdSelector
         {
            typedef Eigen::SparseLU<TMatrix, Eigen::COLAMDOrdering<typename TMatrix::Index> > Type;
         };
      }
   }
#endif //GEOMHDISCC_SPSPDLINALG_SPARSELU

// MUMPS Version for SPD solve
#ifdef GEOMHDISCC_SPSPDLINALG_MUMPS
   // Include the right header
   #include "../External/Interfaces/MumpsLU.hpp"

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers for SPD matries. Sets up the MUMPS version. 
          */
         template<typename TMatrix> struct SparseSpdSelector
         {
            typedef Eigen::MumpsLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPSPDLINALG_MUMPS

// SparseLU Version for SPD solve
#ifdef GEOMHDISCC_SPSPDLINALG_SIMPLICIALLLT
   // Include the right header
   #include <Eigen/SparseCholesky>
   #include <Eigen/OrderingMethods>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers for SPD matries. Sets up Eigen's SimplicialLLT version. 
          */
         template <typename TMatrix> struct SparseSpdSelector
         {
            typedef Eigen::SimplicialLLT<TMatrix, Eigen::Lower, Eigen::AMDOrdering<typename TMatrix::Index> > Type;
         };
      }
   }
#endif //GEOMHDISCC_SPSPDLINALG_SIMPLICIALLLT

// SimplicialLDLT Version for SPD solve
#ifdef GEOMHDISCC_SPSPDLINALG_SIMPLICIALLDLT
   // Include the right header
   #include <Eigen/SparseCholesky>
   #include <Eigen/OrderingMethods>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers for SPD matries. Sets up Eigen's SimplicialLLT version. 
          */
         template <typename TMatrix> struct SparseSpdSelector
         {
            typedef Eigen::SimplicialLDLT<TMatrix, Eigen::Lower, Eigen::AMDOrdering<typename TMatrix::Index> > Type;
         };
      }
   }
#endif //GEOMHDISCC_SPSPDLINALG_SIMPLICIALLDLT

// UmfPack Version for triangular solve
#ifdef GEOMHDISCC_SPTRILINALG_UMFPACK
   // Include the right header
   #include <Eigen/UmfPackSupport>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers for triangular matrices. Sets up the UmfPack version. 
          */
         template <typename TMatrix> struct SparseTriSelector
         {
            typedef Eigen::UmfPackLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPTRILINALG_UMFPACK

// SparseLU Version for triangular solve
#ifdef GEOMHDISCC_SPTRILINALG_SPARSELU
   // Include the right header
   #include <Eigen/SparseLU>
   #include <Eigen/OrderingMethods>

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers for triangular matrices. Sets up Eigen's SparseLU version. 
          */
         template <typename TMatrix> struct SparseTriSelector
         {
            typedef Eigen::SparseLU<TMatrix, Eigen::COLAMDOrdering<typename TMatrix::Index> > Type;
         };
      }
   }
#endif //GEOMHDISCC_SPTRILINALG_SPARSELU

// MUMPS Version for triangular solve
#ifdef GEOMHDISCC_SPTRILINALG_MUMPS
   // Include the right header
   #include "../External/Interfaces/MumpsLU.hpp"

   namespace GeoMHDiSCC {
      namespace Solver {
         /**
          * @brief Selector for the different implementations of the sparse solvers for triangular matrices. Sets up the MUMPS version. 
          */
         template<typename TMatrix> struct SparseTriSelector
         {
            typedef Eigen::MumpsLU<TMatrix> Type;
         };
      }
   }
#endif //GEOMHDISCC_SPTRILINALG_MUMPS

#endif // SPARSESOLVERSELECTOR_HPP
