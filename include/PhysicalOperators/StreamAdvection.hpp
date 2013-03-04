/** \file StreamAdvection.hpp
 *  \brief Implementation of a generic streamfunction advection
 *
 *  \mhdBug Needs test
 */

#ifndef STREAMADVECTION_HPP
#define STREAMADVECTION_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldComponents.hpp"
#include "VectorFields/VectorField.hpp"

namespace GeoMHDiSCC {

namespace Physical {

   /**
    * \brief Implementation of a generic streamfunction advection
    */
   class StreamAdvection
   {
      public:
         /**
          * @brief Set S to streamfunction advection product
          */
         template <int COMPONENTS> static void set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add streamfunction advection product to S
          */
         template <int COMPONENTS> static void add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract streamfunction advection product from S
          */
         template <int COMPONENTS> static void sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         StreamAdvection();

         /**
          * @brief Empty destructor
          */
         ~StreamAdvection();
   };

   template <int COMPONENTS> void StreamAdvection::set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(-c*(dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array()).matrix());
                                                                              
         rS.addData(c*(dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      } else
      {
         rS.setData(-(dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array()).matrix());

         rS.addData((dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      }
   }

   template <int COMPONENTS> void StreamAdvection::add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array()).matrix());

         rS.addData(c*(dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      } else
      {
         rS.subData((dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array()).matrix());

         rS.addData((dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      }
   }

   template <int COMPONENTS> void StreamAdvection::sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array()).matrix());

         rS.subData(c*(dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      } else
      {
         rS.addData((dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array()).matrix());

         rS.subData((dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      }
   }
}
}

#endif // STREAMADVECTION_HPP
