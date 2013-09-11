/** 
 * @file StreamHeatAdvection.hpp
 * @brief Implementation of a generic streamfunction advection including conducting state
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
#include "Enums/FieldIds.hpp"
#include "VectorFields/VectorField.hpp"

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
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         template <int COMPONENTS> static void set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add streamfunction advection product to S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         template <int COMPONENTS> static void add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract streamfunction advection product from S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         template <int COMPONENTS> static void sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);
         
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

   template <int COMPONENTS> void StreamHeatAdvection::set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(-c*(dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array() + dPsi.comp(FieldComponents::Physical::TWO).data().array()).matrix());
                                                                              
         rS.addData(c*(dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      } else
      {
         rS.setData(-(dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array() + dPsi.comp(FieldComponents::Physical::TWO).data().array()).matrix());

         rS.addData((dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      }
   }

   template <int COMPONENTS> void StreamHeatAdvection::add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array() + dPsi.comp(FieldComponents::Physical::TWO).data().array()).matrix());

         rS.addData(c*(dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      } else
      {
         rS.subData((dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array() + dPsi.comp(FieldComponents::Physical::TWO).data().array()).matrix());

         rS.addData((dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      }
   }

   template <int COMPONENTS> void StreamHeatAdvection::sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array() + dPsi.comp(FieldComponents::Physical::TWO).data().array()).matrix());

         rS.subData(c*(dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      } else
      {
         rS.addData((dPsi.comp(FieldComponents::Physical::TWO).data().array()*w.comp(FieldComponents::Physical::ONE).data().array() + dPsi.comp(FieldComponents::Physical::TWO).data().array()).matrix());

         rS.subData((dPsi.comp(FieldComponents::Physical::ONE).data().array()*w.comp(FieldComponents::Physical::TWO).data().array()).matrix());
      }
   }
}
}

#endif // STREAMHEATADVECTION_HPP
