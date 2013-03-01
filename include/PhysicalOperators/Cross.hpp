/** \file Cross.hpp
 *  \brief Implementation of a generic vector cross product
 *
 *  \mhdBug Needs test
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
#include "VectorField/VectorField.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace Physical {

   /**
    * @brief Implementation of a generic vector cross product
    */
   class Cross
   {
      public:
         /**
          * @brief Set S to cross product component
          */
         template <int COMPONENTS> static void set(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add cross product component to S
          */
         template <int COMPONENTS> static void add(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract cross product component from S
          */
         template <int COMPONENTS> static void sub(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);
         
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

   template <int COMPONENTS>  void Cross::set(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      { 
         switch(id)
         {
            case(FieldComponents::ONE):
               rS.setData((v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array() - v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array()).matrix()*c);
               break;
            case(FieldComponents::TWO):
               rS.setData((v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array() - v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array()).matrix()*c);
               break;
            case(FieldComponents::THREE):
               rS.setData((v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array() - v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array()).matrix()*c);
               break;
         }
      } else
      {
         switch(id)
         {
            case(FieldComponents::ONE):
               rS.setData(v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array() - v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array());
               break;
            case(FieldComponents::TWO):
               rS.setData(v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array() - v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array());
               break;
            case(FieldComponents::THREE):
               rS.setData(v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array() - v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array());
               break;
         }
      }
   }

   template <int COMPONENTS> void Cross::add(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         switch(id)
         {
            case(FieldComponents::ONE):
               rS.addData((v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array() - v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array()).matrix()*c);
               break;
            case(FieldComponents::TWO):
               rS.addData((v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array() - v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array()).matrix()*c);
               break;
            case(FieldComponents::THREE):
               rS.addData((v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array() - v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array()).matrix()*c);
               break;
         }
      } else
      {
         switch(id)
         {
            case(FieldComponents::ONE):
               rS.addData(v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array() - v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array());
               break;
            case(FieldComponents::TWO):
               rS.addData(v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array() - v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array());
               break;
            case(FieldComponents::THREE):
               rS.addData(v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array() - v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array());
               break;
         }
      }
   }

   template <int COMPONENTS> void Cross::sub(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         switch(id)
         {
            case(FieldComponents::ONE):
               rS.subData((v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array() - v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array()).matrix()*c);
               break;
            case(FieldComponents::TWO):
               rS.subData((v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array() - v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array()).matrix()*c);
               break;
            case(FieldComponents::THREE):
               rS.subData((v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array() - v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array()).matrix()*c);
               break;
         }
      } else
      {
         switch(id)
         {
            case(FieldComponents::ONE):
               rS.subData(v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array() - v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array());
               break;
            case(FieldComponents::TWO):
               rS.subData(v.comp(FieldsComponents::Physical::THREE).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array() - v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::THREE).data().array());
               break;
            case(FieldComponents::THREE):
               rS.subData(v.comp(FieldsComponents::Physical::ONE).data().array() * w.comp(FieldsComponents::Physical::TWO).data().array() - v.comp(FieldsComponents::Physical::TWO).data().array() * w.comp(FieldsComponents::Physical::ONE).data().array());
               break;
         }
      }
   }
}
}

#endif // CROSS_HPP
