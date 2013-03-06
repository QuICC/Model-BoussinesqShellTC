/** \file FieldTools.hpp
 *  \brief Some useful tools to work with fields 
 */

#ifndef FIELDTOOLS_HPP
#define FIELDTOOLS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <vector>
#include <tr1/tuple>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Single configuration class for the different scalar fields
    */
   class FieldTools
   {
      public:
         /**
          * @brief Create const field data information
          */
         template <typename T, Dimensions::Type DIMENSION, template <typename,Dimensions::Type> class TField> static std::vector<std::tr1::tuple<int, int , const T *> > createInfo(const TField<T,DIMENSION>& field);

         /**
          * @brief Create field data information
          */
         template <typename T, Dimensions::Type DIMENSION, template <typename,Dimensions::Type> class TField> static std::vector<std::tr1::tuple<int, int , T *> >  createInfo(TField<T,DIMENSION>& rField);
         
      protected:

      private:
         /**
          * @brief Constructor 
          */
         FieldTools();

         /**
          * @brief Destructor
          */
         virtual ~FieldTools();
   };

   template <typename T, Dimensions::Type DIMENSION, template <typename,Dimensions::Type> class TField> std::vector<std::tr1::tuple<int, int , const T *> > FieldTools::createInfo(const TField<T,DIMENSION>& field)
   {
      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, const T *> > fieldInfo;

      // Shift for the pointers
      int shift = 0;

      // Loop over all slices
      for(int i = 0; i < field.nSlice(); i++)
      {
         int rows = field.slice(i).rows();
         int cols = field.slice(i).cols();

         fieldInfo.push_back(std::tr1::make_tuple(rows, cols, field.data().data() + shift));

         shift += rows*cols;
      }

      return fieldInfo;
   }

   template <typename T, Dimensions::Type DIMENSION, template <typename,Dimensions::Type> class TField> std::vector<std::tr1::tuple<int, int , T *> > FieldTools::createInfo(TField<T,DIMENSION>& rField)
   {
      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, T *> > fieldInfo;

      // Shift for the pointers
      int shift = 0;

      // Loop over all slices
      for(int i = 0; i < rField.nSlice(); i++)
      {
         int rows = rField.slice(i).rows();
         int cols = rField.slice(i).cols();

         fieldInfo.push_back(std::tr1::make_tuple(rows, cols, rField.rData().data() + shift));

         shift += rows*cols;
      }

      return fieldInfo;
   }

}
}

#endif // FIELDTOOLS_HPP
