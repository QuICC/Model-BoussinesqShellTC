/** 
 * @file FieldTools.hpp
 * @brief Some useful tools to work with fields 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Exceptions/Exception.hpp"
#include "Enums/Dimensions.hpp"
#include "Enums/Arithmetics.hpp"

namespace QuICC {

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

         /**
          * @brief Combine two field with arithmetic operation
          */
         template <typename T, Dimensions::Type DIMENSION, template <typename,Dimensions::Type> class TField> static void combine(TField<T,DIMENSION>& rRhs, const TField<T,DIMENSION>& lhs, const Arithmetics::Id arithId);

         /**
          * @brief Set field to negative
          */
         template <typename T, Dimensions::Type DIMENSION, template <typename,Dimensions::Type> class TField> static void negative(TField<T,DIMENSION>& rRhs);
         
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
      fieldInfo.reserve(field.nSlice());
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
      fieldInfo.reserve(rField.nSlice());
      for(int i = 0; i < rField.nSlice(); i++)
      {
         int rows = rField.slice(i).rows();
         int cols = rField.slice(i).cols();

         fieldInfo.push_back(std::tr1::make_tuple(rows, cols, rField.rData().data() + shift));

         shift += rows*cols;
      }

      return fieldInfo;
   }

   template <typename T, Dimensions::Type DIMENSION, template <typename,Dimensions::Type> class TField> void FieldTools::combine(TField<T,DIMENSION>& rRhs, const TField<T,DIMENSION>& lhs, const Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::ADD || arithId == Arithmetics::SUB);

      if(arithId == Arithmetics::ADD)
      {
         rRhs.addData(lhs.data());
      } else if(arithId == Arithmetics::SUB)
      {
         rRhs.subData(lhs.data());
      } else
      {
         throw Exception("Unknown arithmetic operation in combine!");
      }
   }

   template <typename T, Dimensions::Type DIMENSION, template <typename,Dimensions::Type> class TField> void FieldTools::negative(TField<T,DIMENSION>& rRhs)
   {
      rRhs.rData() = -rRhs.data();
   }

}
}

#endif // FIELDTOOLS_HPP
