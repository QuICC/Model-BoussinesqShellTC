/**
 * @file MpiTypes.hpp
 * @brief Definition of a simple MPI datatype converters
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// Only define in MPI case
#ifdef GEOMHDISCC_MPI

#ifndef MPITYPES_HPP
#define MPITYPES_HPP

// System includes
//
#include <complex>
#include <mpi.h>

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

// Only define in MPI case
#endif //GEOMHDISCC_MPI
