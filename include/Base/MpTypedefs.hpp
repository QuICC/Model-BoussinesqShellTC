/**
 * @file MpTypedefs.hpp
 * @brief Definition of typedefs for muliple precision computations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MPTYPEDEFS_HPP
#define MPTYPEDEFS_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//
#include <Eigen/Core>
#include <Eigen/MPRealSupport>

// Project includes
//

namespace GeoMHDiSCC {

   /**
    * @name Basic scalar types typedefs
    */
   //@{
   /// Typedef for multiple precision floating point type value
   typedef mpfr::mpreal MHDMpFloat;
   //@}

   /**
    * @name Array types typedefs
    */
   //@{
   /// Typedef for an array of multiple precision float values
   typedef Eigen::Matrix<MHDMpFloat, Eigen::Dynamic, 1>   MpArray;
   //@}

   /**
    * @name Matrix types typedefs
    */
   //@{
   /// Typedef for a matrix of multiple precision float values
   typedef Eigen::Matrix<MHDMpFloat, Eigen::Dynamic, Eigen::Dynamic>  MpMatrix;
   //@}

   /**
    * @name Shared pointer typedefs
    */
   //@{
   /// Typedef for an smart reference counting pointer on an array of multiple precision real values
   typedef SharedPtrMacro<MpArray>   SharedMpArray;
   /// Typedef for an smart reference counting pointer on an matrix of real values
   typedef SharedPtrMacro<MpMatrix>   SharedMpMatrix;
   //@}

}

#endif // MPTYPEDEFS_HPP
