/** \file MpiTypes.hpp
 *  \brief Implementation of a simple MPI datatype converter
 */

#ifndef MPITYPES_HPP
#define MPITYPES_HPP

// System includes
//
#include <complex>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief This class implements a few type conversion routines for the MPI library
    */
   struct MpiTypes
   {
      /**
       * @brief Get MPI predefined type
       */
      template <typename T> static MPI_Datatype type();
   };

   /**
    * @brief Specialised method for integer type
    */
   template <> inline MPI_Datatype MpiTypes::type<int>()
   {
      return MPI_INT;
   }

   /**
    * @brief Specialised method for float type
    */
   template <> inline MPI_Datatype MpiTypes::type<float>()
   {
      return MPI_FLOAT;
   }

   /**
    * @brief Specialised method for double type
    */
   template <> inline MPI_Datatype MpiTypes::type<double>()
   {
      return MPI_DOUBLE;
   }

   /**
    * @brief Specialised method for a complex<double>
    */
   template <> inline MPI_Datatype MpiTypes::type<std::complex<double> >()
   {
      return MPI_DOUBLE_COMPLEX;
   }

}
}

#endif // MPITYPES_HPP
