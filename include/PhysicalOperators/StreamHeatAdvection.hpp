/** \file StreamHeatAdvection.hpp
 *  \brief Implementation of a generic streamfunction advection including conducting state
 *
 *  \mhdBug Needs test
 */

#ifndef STREAMHEATADVECTION_HPP
#define STREAMHEATADVECTION_HPP

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

namespace GeoMHDiSCC {

namespace Physical {

   /**
    * \brief Implementation of a generic streamfunction advection including conducting state
    */
   class StreamHeatAdvection
   {
      public:
         /**
          * @brief Set S to streamfunction advection product
          */
         template <int COMPONENTS> static void set(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add streamfunction advection product to S
          */
         template <int COMPONENTS> static void add(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract streamfunction advection product from S
          */
         template <int COMPONENTS> static void sub(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         StreamHeatAdvection();

         /**
          * @brief Empty destructor
          */
         ~StreamHeatAdvection();
   };

   template <int COMPONENTS> void StreamHeatAdvection::set(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(-c*(dPsi.comp(FieldsComponents::Physical::TWO).data().array()*w.comp(FieldsComponents::Physical::ONE).data().array() + dPsi.comp(FieldsComponents::Physical::TWO).data().array()).matrix());
                                                                              
         rS.addData(c*(dPsi.comp(FieldsComponents::Physical::ONE).data().array()*w.comp(FieldsComponents::Physical::TWO).data().array()).matrix());
      } else
      {
         rS.setData(-(dPsi.comp(FieldsComponents::Physical::TWO).data().array()*w.comp(FieldsComponents::Physical::ONE).data().array() + dPsi.comp(FieldsComponents::Physical::TWO).data().array()).matrix());

         rS.addData((dPsi.comp(FieldsComponents::Physical::ONE).data().array()*w.comp(FieldsComponents::Physical::TWO).data().array()).matrix());
      }
   }

   template <int COMPONENTS> void StreamHeatAdvection::add(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(dPsi.comp(FieldsComponents::Physical::TWO).data().array()*w.comp(FieldsComponents::Physical::ONE).data().array() + dPsi.comp(FieldsComponents::Physical::TWO).data().array()).matrix());

         rS.addData(c*(dPsi.comp(FieldsComponents::Physical::ONE).data().array()*w.comp(FieldsComponents::Physical::TWO).data().array()).matrix());
      } else
      {
         rS.subData((dPsi.comp(FieldsComponents::Physical::TWO).data().array()*w.comp(FieldsComponents::Physical::ONE).data().array() + dPsi.comp(FieldsComponents::Physical::TWO).data().array()).matrix());

         rS.addData((dPsi.comp(FieldsComponents::Physical::ONE).data().array()*w.comp(FieldsComponents::Physical::TWO).data().array()).matrix());
      }
   }

   template <int COMPONENTS> void StreamHeatAdvection::sub(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(dPsi.comp(FieldsComponents::Physical::TWO).data().array()*w.comp(FieldsComponents::Physical::ONE).data().array() + dPsi.comp(FieldsComponents::Physical::TWO).data().array()).matrix());

         rS.subData(c*(dPsi.comp(FieldsComponents::Physical::ONE).data().array()*w.comp(FieldsComponents::Physical::TWO).data().array()).matrix());
      } else
      {
         rS.addData((dPsi.comp(FieldsComponents::Physical::TWO).data().array()*w.comp(FieldsComponents::Physical::ONE).data().array() + dPsi.comp(FieldsComponents::Physical::TWO).data().array()).matrix());

         rS.subData((dPsi.comp(FieldsComponents::Physical::ONE).data().array()*w.comp(FieldsComponents::Physical::TWO).data().array()).matrix());
      }
   }
}
}

#endif // STREAMHEATADVECTION_HPP
