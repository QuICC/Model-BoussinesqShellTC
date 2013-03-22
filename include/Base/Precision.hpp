/** \file Precision.hpp
 *  \brief Definition of small class containing normal or MP typedefs for internal computations
 */

#ifndef PRECISION_HPP
#define PRECISION_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

#ifdef GEOMHDISCC_MULTPRECISION
   #include "Base/MpTypedefs.hpp"
   
   /// Define a small macro to replace float constants to strings in the case of MP computations
   #define MHD_MP(c) #c
#else
   /// For normal computations the macro does nothing
   #define MHD_MP(c) c
#endif // GEOMHDISCC_MULTPRECISION

namespace GeoMHDiSCC {

   namespace internal {

      #ifdef GEOMHDISCC_MULTPRECISION
      //////////////////////////////////////////////////////////////////         
         /// Typedef for the internal float type
         typedef GeoMHDiSCC::MHDMpFloat MHDFloat;

         /// Typedef for the internal Array type
         typedef GeoMHDiSCC::MpArray Array;

         /// Typedef for the internal Matrix type
         typedef GeoMHDiSCC::MpMatrix Matrix;

         /// Typedef for the smart internal Array type
         typedef GeoMHDiSCC::SharedMpArray SharedArray;

         /// Typedef for the smart internal Matrix type
         typedef GeoMHDiSCC::SharedMpMatrix SharedMatrix;
      //////////////////////////////////////////////////////////////////
      #else
      //////////////////////////////////////////////////////////////////
         /// Typedef for the internal float type
         typedef GeoMHDiSCC::MHDFloat MHDFloat;

         /// Typedef for the internal Array type
         typedef GeoMHDiSCC::Array Array;

         /// Typedef for the internal Matrix type
         typedef GeoMHDiSCC::Matrix Matrix;

         /// Typedef for the smart internal Array type
         typedef GeoMHDiSCC::SharedArray SharedArray;

         /// Typedef for the smart internal Matrix type
         typedef GeoMHDiSCC::SharedMatrix SharedMatrix;
      //////////////////////////////////////////////////////////////////
      #endif // GEOMHDISCC_MULTPRECISION
   }

   /**
    * @brief Simple class holding some typedefs to allow for internal MP computations
    */
   class Precision
   {
      public:

         /**
          * @brief Precision dependent mathematical constant \f$\pi\f$
          */
         static const internal::MHDFloat PI;

         /**
          * @brief Initialise the precision setup
          */
         static void init();

         /**
          * @brief Cast the internal smart Array to an external one
          *
          * @param spIArr Internal smart Array to cast
          */
         static SharedArray cast(internal::SharedArray spIArr);

         /**
          * @brief Cast the internal smart Matrix to an external one
          *
          * @param spIMat Internal smart Matrix to cast
          */
         static SharedMatrix cast(internal::SharedMatrix spIMat);

         /**
          * @brief Cast the internal smart Array to an external one
          *
          * @param rIArr Internal Array to cast
          */
         static Array cast(const internal::Array& rIArr);

         /**
          * @brief Cast the internal smart Matrix to an external one
          *
          * @param rIMat Internal Matrix to cast
          */
         static Matrix cast(const internal::Matrix& rIMat);

      private:
         /**
         * @brief Empty constructor
         */
         Precision();

         /**
         * @brief Simple empty destructor
         */
         ~Precision();
   };

   inline void Precision::init()
   {
      #ifdef GEOMHDISCC_MULTPRECISION
         mpfr::mpreal::set_default_prec(256);
      #endif // GEOMHDISCC_MULTPRECISION
   }

   inline SharedArray Precision::cast(internal::SharedArray spIArr)
   {
      #ifdef GEOMHDISCC_MULTPRECISION
         SharedArray spArr(new Array(spIArr->size()));

         // Loop over whole array
         for(int i=0; i < spIArr->size(); i++)
         {
            (*spArr)(i) = (*spIArr)(i).toDouble();
         }

         return spArr;
      #else
         return spIArr;
      #endif // GEOMHDISCC_MULTPRECISION
   }

   inline SharedMatrix Precision::cast(internal::SharedMatrix spIMat)
   {
      #ifdef GEOMHDISCC_MULTPRECISION
         SharedMatrix spMat(new Matrix(spIMat->rows(),spIMat->cols()));

         // Loop over whole matrix
         for(int j=0; j < spIMat->cols(); j++)
         {
            for(int i=0; i < spIMat->rows(); i++)
            {
               (*spMat)(i,j) = (*spIMat)(i,j).toDouble();
            }
         }

         return spMat;
      #else
         return spIMat;
      #endif // GEOMHDISCC_MULTPRECISION
   }

   inline Array Precision::cast(const internal::Array& rIArr)
   {
      #ifdef GEOMHDISCC_MULTPRECISION
         Array arr(rIArr.size());

         // Loop over whole array
         for(int i=0; i < rIArr.size(); i++)
         {
            arr(i) = rIArr(i).toDouble();
         }

         return arr;
      #else
         return rIArr;
      #endif // GEOMHDISCC_MULTPRECISION
   }

   inline Matrix Precision::cast(const internal::Matrix& rIMat)
   {
      #ifdef GEOMHDISCC_MULTPRECISION
         Matrix mat(rIMat.rows(),rIMat.cols());

         // Loop over whole matrix
         for(int j=0; j < rIMat.cols(); j++)
         {
            for(int i=0; i < rIMat.rows(); i++)
            {
               mat(i,j) = rIMat(i,j).toDouble();
            }
         }

         return mat;
      #else
         return rIMat;
      #endif // GEOMHDISCC_MULTPRECISION
   }

#ifdef GEOMHDISCC_MULTPRECISION
   /// Create a namespace alias for the internal precision stuff pointing to mpfr namespace
   namespace  precision = mpfr;
#else
   /// Create a namespace alias for the internal precision stuff pointing to std namespace
   namespace  precision = std;
#endif // GEOMHDISCC_MULTPRECISION
}

#endif // PRECISION_HPP
