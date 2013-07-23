/** \file SparseSolverMacro.h
 *  \brief Preprocessor macros used to select different sparse solvers.
 */

#ifndef SPARSESOLVERMACRO_H
#define SPARSESOLVERMACRO_H

// SuperLU sparse solver
#ifdef GEOMHDISCC_SPLINALG_SUPERLU
   // Include the right header
   #include <Eigen/SuperLUSupport>

   /**
    * @def SparseSolverMacro
    * Macro allowing to use different implementations of the sparse solvers.
    * Sets up the SuperLU version.
    */
   #define SparseSolverMacro Eigen::SuperLU
#endif //GEOMHDISCC_SPLINALG_SUPERLU

// UMFPack Version
#ifdef GEOMHDISCC_SPLINALG_UMFPACK
   // Include the right header
   #include <Eigen/UmfPackSupport>

   /**
    * @def SparseSolverMacro
    * Macro allowing to use different implementations of the sparse solvers.
    * Sets up the UMFPack version.
    */
   #define SparseSolverMacro Eigen::UmfPackLU
#endif //GEOMHDISCC_SPLINALG_UMFPACK

// SparseLU Version
#ifdef GEOMHDISCC_SPLINALG_SPARSELU
   // Include the right header
   #include <Eigen/SparseLU>

   /**
    * @def SparseSolverMacro
    * Macro allowing to use different implementations of the sparse solvers.
    * Sets up the Eigen internal SparseLU version.
    */
   #define SparseSolverMacro Eigen::SparseLU
#endif //GEOMHDISCC_SPLINALG_SPARSELU

#endif // SPARSESOLVERMACRO_H
