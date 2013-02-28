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
#include "VectorField/VectorField.hpp"

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
         StreamAdvection();

         /**
          * @brief Empty destructor
          */
         ~StreamAdvection();
   };

   template <int COMPONENTS> void StreamAdvection::set(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.rData() = -c*(dPsi.comp(1).data().array()*w.comp(0).data().array()).matrix();
                                                                              
         rS.rData() += c*(dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      } else
      {
         rS.rData() = -(dPsi.comp(1).data().array()*w.comp(0).data().array()).matrix();

         rS.rData() += (dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      }
   }

   template <int COMPONENTS> void StreamAdvection::add(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.rData() -= c*(dPsi.comp(1).data().array()*w.comp(0).data().array()).matrix();

         rS.rData() += c*(dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      } else
      {
         rS.rData() -= (dPsi.comp(1).data().array()*w.comp(0).data().array()).matrix();

         rS.rData() += (dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      }
   }

   template <int COMPONENTS> void StreamAdvection::sub(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.rData() += c*(dPsi.comp(1).data().array()*w.comp(0).data().array()).matrix();

         rS.rData() -= c*(dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      } else
      {
         rS.rData() += (dPsi.comp(1).data().array()*w.comp(0).data().array()).matrix();

         rS.rData() -= (dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      }
   }
}
}

#endif // STREAMADVECTION_HPP
