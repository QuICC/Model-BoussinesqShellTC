/** 
 * @file StreamAdvection.hpp
 * @brief Implementation of a generic streamfunction advection
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Enums/FieldIds.hpp"
#include "VectorFields/VectorField.hpp"

namespace GeoMHDiSCC {

namespace Physical {

   /**
    * @brief Implementation of a generic streamfunction advection
    */
   class StreamAdvection
   {
      public:
         /**
          * @brief Set S to streamfunction advection product
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)q = -\partial_y\psi\partial_x q + \partial_x\psi\partial_y q\f$
          */
         template <int COMPONENTS> static void set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0, FieldComponents::Physical::Id compX = FieldComponents::Physical::ONE, FieldComponents::Physical::Id compY = FieldComponents::Physical::TWO);

         /**
          * @brief Add streamfunction advection product to S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)q = -\partial_y\psi\partial_x q + \partial_x\psi\partial_y q\f$
          */
         template <int COMPONENTS> static void add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0, FieldComponents::Physical::Id compX = FieldComponents::Physical::ONE, FieldComponents::Physical::Id compY = FieldComponents::Physical::TWO);

         /**
          * @brief Substract streamfunction advection product from S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)q = -\partial_y\psi\partial_x q + \partial_x\psi\partial_y q\f$
          */
         template <int COMPONENTS> static void sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0, FieldComponents::Physical::Id compX = FieldComponents::Physical::ONE, FieldComponents::Physical::Id compY = FieldComponents::Physical::TWO);
         
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

   template <int COMPONENTS> void StreamAdvection::set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c, FieldComponents::Physical::Id compX, FieldComponents::Physical::Id compY)
   {
      if(c != 1.0)
      {
         rS.setData(-c*(dPsi.comp(compY).data().array()*w.comp(compX).data().array()).matrix());
                                                                              
         rS.addData(c*(dPsi.comp(compX).data().array()*w.comp(compY).data().array()).matrix());
      } else
      {
         rS.setData(-(dPsi.comp(compY).data().array()*w.comp(compX).data().array()).matrix());

         rS.addData((dPsi.comp(compX).data().array()*w.comp(compY).data().array()).matrix());
      }
   }

   template <int COMPONENTS> void StreamAdvection::add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c, FieldComponents::Physical::Id compX, FieldComponents::Physical::Id compY)
   {
      if(c != 1.0)
      {
         rS.subData(c*(dPsi.comp(compY).data().array()*w.comp(compX).data().array()).matrix());

         rS.addData(c*(dPsi.comp(compX).data().array()*w.comp(compY).data().array()).matrix());
      } else
      {
         rS.subData((dPsi.comp(compY).data().array()*w.comp(compX).data().array()).matrix());

         rS.addData((dPsi.comp(compX).data().array()*w.comp(compY).data().array()).matrix());
      }
   }

   template <int COMPONENTS> void StreamAdvection::sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c, FieldComponents::Physical::Id compX, FieldComponents::Physical::Id compY)
   {
      if(c != 1.0)
      {
         rS.addData(c*(dPsi.comp(compY).data().array()*w.comp(compX).data().array()).matrix());

         rS.subData(c*(dPsi.comp(compX).data().array()*w.comp(compY).data().array()).matrix());
      } else
      {
         rS.addData((dPsi.comp(compY).data().array()*w.comp(compX).data().array()).matrix());

         rS.subData((dPsi.comp(compX).data().array()*w.comp(compY).data().array()).matrix());
      }
   }
}
}

#endif // STREAMADVECTION_HPP
