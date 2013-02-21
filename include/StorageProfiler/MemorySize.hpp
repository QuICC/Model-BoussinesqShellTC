/** \file MemorySize.hpp
 *  \brief Simple template for the memory usage of the different types
 */

#ifndef MEMORYSIZE_HPP
#define MEMORYSIZE_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Debug {

   /**
    * @brief Simple template for the memory usage of the different types
    */
   template <typename T> struct MemorySize
   {
   };

   template <> struct MemorySize<char>
   {
      static const int BYTES = 1;
   };

   template <> struct MemorySize<int>
   {
      static const int BYTES = 4;
   };

   template <> struct MemorySize<float>
   {
      static const int BYTES = 4;
   };

   template <> struct MemorySize<double>
   {
      static const int BYTES = 8;
   };

   template <> struct MemorySize<std::complex<float> >
   {
      static const int BYTES = 8;
   };

   template <> struct MemorySize<std::complex<double> >
   {
      static const int BYTES = 16;
   };

}
}

#endif // MEMORYSIZE_HPP
