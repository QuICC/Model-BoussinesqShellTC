/** \file Dot.hpp
 *  \brief Implementation of a generic scalar product
 */

#ifndef DOT_HPP
#define DOT_HPP

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

namespace GeoMHDiSCC {

namespace Physical {

   /**
    * @brief Implementation of a generic scalar product
    */
   class Dot
   {
      public:
         /**
          * @brief Set S to scalar product
          */
         template <int COMPONENTS> static void set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add scalar product to S
          */
         template <int COMPONENTS> static void add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract scalar product from S
          */
         template <int COMPONENTS> static void sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         Dot();

         /**
          * @brief Empty destructor
          */
         ~Dot();
   };

   template <int COMPONENTS> void Dot::set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(v.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::ONE).data().array()).matrix());

         for(int i = 1; i < COMPONENTS; i++)
         {
            rS.addData(c*(v.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()*w.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()).matrix());
         }
      } else
      {
         rS.setData((v.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::ONE).data().array()).matrix());

         for(int i = 1; i < COMPONENTS; i++)
         {
            rS.addData((v.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()*w.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()).matrix());
         }
      }
   }

   template <int COMPONENTS> void Dot::add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         for(int i = 0; i < COMPONENTS; i++)
         {
            rS.addData(c*(v.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()*w.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()).matrix());
         }
      } else
      {
         for(int i = 0; i < COMPONENTS; i++)
         {
            rS.addData((v.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()*w.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()).matrix());
         }
      }
   }

   template <int COMPONENTS> void Dot::sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         for(int i = 0; i < COMPONENTS; i++)
         {
            rS.subData(c*(v.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()*w.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()).matrix());
         }
      } else
      {
         for(int i = 0; i < COMPONENTS; i++)
         {
            rS.subData((v.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()*w.comp(static_cast<FieldComponents::Physical::Id>(i)).data().array()).matrix());
         }
      }
   }
}
}

#endif // DOT_HPP
