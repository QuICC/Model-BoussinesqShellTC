/** 
 * @file StreamHeatAdvection.hpp
 * @brief Implementation of a generic streamfunction advection including conducting state
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of a generic streamfunction advection including conducting state
    */
   template <FieldComponents::Physical::Id TXComp = FieldComponents::Physical::ONE, FieldComponents::Physical::Id TYComp = FieldComponents::Physical::TWO> class StreamHeatAdvection
   {
      public:
         /**
          * @brief Set S to streamfunction advection product
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         static void set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add streamfunction advection product to S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         static void add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract streamfunction advection product from S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         static void sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);
         
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

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> void StreamHeatAdvection<TXComp,TYComp>::set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());
                                                                              
         rS.subData(c*(dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      } else
      {
         rS.setData((dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());
                                                                              
         rS.subData((dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> void StreamHeatAdvection<TXComp,TYComp>::add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData(c*(dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      } else
      {
         rS.addData((dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData((dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> void StreamHeatAdvection<TXComp,TYComp>::sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.addData(c*(dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      } else
      {
         rS.subData((dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.addData((dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      }
   }
}
}

#endif // STREAMHEATADVECTION_HPP
