/** 
 * @file VelocityHeatAdvection.hpp
 * @brief Implementation of a generic primitive velocity heat advection:
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
    * @brief Implementation of a generic primitive velocity heat advection in 3D space
    */
   template <FieldComponents::Physical::Id TONE = FieldComponents::Physical::ONE, FieldComponents::Physical::Id TTWO = FieldComponents::Physical::TWO, FieldComponents::Physical::Id TTHREE = FieldComponents::Physical::THREE> class VelocityHeatAdvection
   {
      public:
         /**
          * @brief Set S to primitive velocity heat advection product
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q - u_z = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q - u_z\f$
          */
          static void set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add primitive velocity advection heat product to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q - u_z\f$
          */
          static void add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract primitive velocity heat advection product from S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q - u_z\f$
          */
          static void sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);
         
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

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> void VelocityHeatAdvection<TONE,TTWO,TTHREE>::set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
                                                                              
         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array() - u.comp(TTHREE).data().array()).matrix());
      } else
      {
         rS.setData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array() - u.comp(TTHREE).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> void VelocityHeatAdvection<TONE,TTWO,TTHREE>::add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
                                                                              
         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array() - u.comp(TTHREE).data().array()).matrix());
      } else
      {
         rS.addData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array() - u.comp(TTHREE).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> void VelocityHeatAdvection<TONE,TTWO,TTHREE>::sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
                                                                              
         rS.subData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.subData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array() - u.comp(TTHREE).data().array()).matrix());
      } else
      {
         rS.subData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.subData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.subData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array() - u.comp(TTHREE).data().array()).matrix());
      }
   }

   /**
    * @brief Implementation of a generic primitive velocity heat advection in 2D space
    */
   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO> class VelocityHeatAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>
   {
      public:
         /**
          * @brief Set S to primitive velocity heat advection product
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          static void set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add primitive velocity advection heat product to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          static void add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract primitive velocity heat advection product from S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          static void sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);
         
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

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO> void VelocityHeatAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::set(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
                                                                              
         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array() - u.comp(TTWO).data().array()).matrix());
      } else
      {
         rS.setData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array() - u.comp(TTWO).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO> void VelocityHeatAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::add(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
                                                                              
         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array() - u.comp(TTWO).data().array()).matrix());
      } else
      {
         rS.addData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array() - u.comp(TTWO).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO> void VelocityHeatAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::sub(Datatypes::PhysicalScalarType &rS, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
                                                                              
         rS.subData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array() - u.comp(TTWO).data().array()).matrix());
      } else
      {
         rS.subData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.subData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array() - u.comp(TTWO).data().array()).matrix());
      }
   }
}
}

#endif // VELOCITYHEATADVECTION_HPP
