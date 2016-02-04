/** 
 * @file Hdf5Types.hpp
 * @brief Implementation of a simple HDF5 datatype converter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef HDF5TYPES_HPP
#define HDF5TYPES_HPP

// System includes
//
#include <complex>
#include <hdf5.h>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace IoHdf5 {

   /**
    * @brief This class implements a few type conversion routines for the HDF5 library
    */
   struct Hdf5Types
   {
      /**
       * @brief Get HDF5 predefined type for a scalar type
       *
       * \tparam T type to convert to HDF5 type
       */
      template <typename T> static hid_t type();
   };

   /**
    * @brief Specialised method for integer type
    */
   template <> inline hid_t Hdf5Types::type<int>()
   {
      return H5T_NATIVE_INT;
   }

   /**
    * @brief Specialised method for float type
    */
   template <> inline hid_t Hdf5Types::type<float>()
   {
      return H5T_NATIVE_FLOAT;
   }

   /**
    * @brief Specialised method for double type
    */
   template <> inline hid_t Hdf5Types::type<double>()
   {
      return H5T_NATIVE_DOUBLE;
   }

   /**
    * @brief Specialised method for a complex<float>
    */
   template <> inline hid_t Hdf5Types::type<std::complex<float> >()
   {
      hsize_t dims = 2;
      return H5Tarray_create(H5T_NATIVE_FLOAT, 1, &dims);
   }

   /**
    * @brief Specialised method for a complex<double>
    */
   template <> inline hid_t Hdf5Types::type<std::complex<double> >()
   {
//      hsize_t dims = 2;
//      return H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &dims);

      std::complex<double> tmp(0,0);
      double d = 0;
      hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof tmp);
      H5Tinsert (complex_id, "r", 0, H5T_NATIVE_DOUBLE);
      H5Tinsert (complex_id, "i", sizeof d, H5T_NATIVE_DOUBLE);

      return complex_id;
   }

}
}

#endif // HDF5TYPES_HPP
