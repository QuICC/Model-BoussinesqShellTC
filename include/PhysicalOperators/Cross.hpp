/** 
 * @file Cross.hpp
 * @brief Implementation of a generic vector cross product
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CROSS_HPP
#define CROSS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "VectorFields/VectorField.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace Physical {

   /**
    * @brief Implementation of a generic vector cross product
    */
   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> class Cross
   {
      public:
         /**
          * @brief Set S to cross product component
          */
         static void set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add cross product component to S
          */
         static void add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract cross product component from S
          */
         static void sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         Cross();

         /**
          * @brief Empty destructor
          */
         ~Cross();
   };

   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> inline void Cross<TFIRST,TSECOND>::set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      { 
         rS.setData(c*(v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.subData(c*(v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      } else
      {
         rS.setData((v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.subData((v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> inline void Cross<TFIRST,TSECOND>::add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.subData(c*(v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      } else
      {
         rS.addData((v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.subData((v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> inline void Cross<TFIRST,TSECOND>::sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.addData(c*(v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      } else
      {
         rS.subData((v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.addData((v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      }
   }
}
}

#endif // CROSS_HPP
