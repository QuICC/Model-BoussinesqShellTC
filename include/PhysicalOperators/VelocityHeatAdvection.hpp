/** 
 * @file VelocityHeatAdvection.hpp
 * @brief Implementation of a generic primitive velocity heat advection
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VELOCITYHEATADVECTION_HPP
#define VELOCITYHEATADVECTION_HPP

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
    * @brief Implementation of a generic primitive velocity heat advection
    */
   template <FieldComponents::Physical::Id TXComp = FieldComponents::Physical::ONE, FieldComponents::Physical::Id TYComp = FieldComponents::Physical::TWO, FieldComponents::Physical::Id TZComp = FieldComponents::Physical::THREE> class VelocityHeatAdvection
   {
      public:
         /**
          * @brief Set S to primitive velocity heat advection product
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          template <int COMPONENTS> static void set(Datatypes::PhysicalScalarType &rS, const Datatypes::PhysicalScalarType &uX, const Datatypes::PhysicalScalarType &uY, const Datatypes::PhysicalScalarType &uZ, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add primitive velocity advection heat product to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          template <int COMPONENTS> static void add(Datatypes::PhysicalScalarType &rS, const Datatypes::PhysicalScalarType &uX, const Datatypes::PhysicalScalarType &uY, const Datatypes::PhysicalScalarType &uZ, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract primitive velocity heat advection product from S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          template <int COMPONENTS> static void sub(Datatypes::PhysicalScalarType &rS, const Datatypes::PhysicalScalarType &uX, const Datatypes::PhysicalScalarType &uY, const Datatypes::PhysicalScalarType &uZ, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         VelocityHeatAdvection();

         /**
          * @brief Empty destructor
          */
         ~VelocityHeatAdvection();
   };

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp, FieldComponents::Physical::Id TZComp> template <int COMPONENTS> void VelocityHeatAdvection<TXComp,TYComp,TZComp>::set(Datatypes::PhysicalScalarType &rS, const Datatypes::PhysicalScalarType &uX, const Datatypes::PhysicalScalarType &uY, const Datatypes::PhysicalScalarType &uZ, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(uX.data().array()*gradQ.comp(TXComp).data().array()).matrix());
                                                                              
         rS.addData(c*(uY.data().array()*gradQ.comp(TYComp).data().array()).matrix());

         rS.addData(c*(uZ.data().array()*gradQ.comp(TZComp).data().array() + uZ.data().array()).matrix());
      } else
      {
         rS.setData((uX.data().array()*gradQ.comp(TXComp).data().array()).matrix());

         rS.addData((uY.data().array()*gradQ.comp(TYComp).data().array()).matrix());

         rS.addData((uZ.data().array()*gradQ.comp(TZComp).data().array() + uZ.data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp, FieldComponents::Physical::Id TZComp> template <int COMPONENTS> void VelocityHeatAdvection<TXComp,TYComp,TZComp>::add(Datatypes::PhysicalScalarType &rS, const Datatypes::PhysicalScalarType &uX, const Datatypes::PhysicalScalarType &uY, const Datatypes::PhysicalScalarType &uZ, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(uX.data().array()*gradQ.comp(TXComp).data().array()).matrix());
                                                                              
         rS.addData(c*(uY.data().array()*gradQ.comp(TYComp).data().array()).matrix());

         rS.addData(c*(uZ.data().array()*gradQ.comp(TZComp).data().array() + uZ.data().array()).matrix());
      } else
      {
         rS.addData((uX.data().array()*gradQ.comp(TXComp).data().array()).matrix());

         rS.addData((uY.data().array()*gradQ.comp(TYComp).data().array()).matrix());

         rS.addData((uZ.data().array()*gradQ.comp(TZComp).data().array() + uZ.data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp, FieldComponents::Physical::Id TZComp> template <int COMPONENTS> void VelocityHeatAdvection<TXComp,TYComp,TZComp>::sub(Datatypes::PhysicalScalarType &rS, const Datatypes::PhysicalScalarType &uX, const Datatypes::PhysicalScalarType &uY, const Datatypes::PhysicalScalarType &uZ, const Datatypes::VectorField<Datatypes::PhysicalScalarType, COMPONENTS, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(uX.data().array()*gradQ.comp(TXComp).data().array()).matrix());
                                                                              
         rS.subData(c*(uY.data().array()*gradQ.comp(TYComp).data().array()).matrix());

         rS.subData(c*(uZ.data().array()*gradQ.comp(TZComp).data().array()+ uZ.data().array()).matrix());
      } else
      {
         rS.subData((uX.data().array()*gradQ.comp(TXComp).data().array()).matrix());

         rS.subData((uY.data().array()*gradQ.comp(TYComp).data().array()).matrix());

         rS.subData((uZ.data().array()*gradQ.comp(TZComp).data().array() + uZ.data().array()).matrix());
      }
   }
}
}

#endif // VELOCITYHEATADVECTION_HPP
