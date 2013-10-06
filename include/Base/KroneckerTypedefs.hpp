/**
 * @file KroneckerTypedefs.hpp
 * @brief Some general typedefs used for kronecker products
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef KRONECKERTYPEDEFS_HPP
#define KRONECKERTYPEDEFS_HPP

// Configuration includes
//

// System includes
//
#include <tr1/tuple>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

   /**
    * @name Kronecker product of real operators
    */
   //@{
   // No eigen dimension
   typedef std::tr1::tuple<SparseMatrix,SparseMatrix,SparseMatrix> KronNoEigenRProduct;
   // Single eigen dimension
   typedef std::tr1::tuple<SparseMatrix,SparseMatrix> KronEigen1DRProduct;
   // Two eigen dimension
   typedef SparseMatrix KronEigen2DRProduct;
   // Three eigen dimension
   typedef MHDFloat KronEigen3DRProduct;
   //@}

   /**
    * @name Sum of Kronecker products of real operators
    */
   //@{
   // No eigen dimension
   typedef std::vector<KronNoEigenRProduct> KronNoEigenRSum;
   // Single eigen dimension
   typedef std::vector<KronEigen1DRProduct> KronEigen1DRSum;
   // Two eigen dimension
   typedef SparseMatrix KronEigen2DRSum;
   // Three eigen dimension
   typedef MHDFloat KronEigen3DRSum;
   //@}

   /**
    * @name Kronecker product of complex operators in decoupled form
    */
   //@{
   // No eigen dimension
   typedef std::tr1::tuple<DecoupledZSparse,DecoupledZSparse,DecoupledZSparse> KronNoEigenZProduct;
   // Single eigen dimension
   typedef std::tr1::tuple<DecoupledZSparse,DecoupledZSparse> KronEigen1DZProduct;
   // Two eigen dimension
   typedef DecoupledZSparse KronEigen2DZProduct;
   // Three eigen dimension
   typedef MHDComplex KronEigen3DZProduct;
   //@}

   /**
    * @name Sum of Kronecker products of complex operators in decoupled form
    */
   //@{
   // No eigen dimension
   typedef std::vector<KronNoEigenZProduct> KronNoEigenZSum;
   // Single eigen dimension
   typedef std::vector<KronEigen1DZProduct> KronEigen1DZSum;
   // Two eigen dimension
   typedef DecoupledZSparse KronEigen2DZSum;
   // Three eigen dimension
   typedef MHDComplex KronEigen3DZSum;
   //@}
}

#endif // KRONECKERTYPEDEFS_HPP