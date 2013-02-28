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
         template <int COMPONENTS> static void set(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add cross product component to S
          */
         template <int COMPONENTS> static void add(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract cross product component from S
          */
         template <int COMPONENTS> static void sub(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c = 1.0);
         
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

   template <int COMPONENTS>  void Cross::set(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c)
   {
      if(c != 1.0)
      { 
         switch(static_cast<int>(id))
         {
            case(0):
               rS.rData() = (v.comp(1).data().array() * w.comp(2).data().array() - v.comp(2).data().array() * w.comp(1).data().array()).matrix()*c;
               break;
            case(1):
               rS.rData() = (v.comp(2).data().array() * w.comp(0).data().array() - v.comp(0).data().array() * w.comp(2).data().array()).matrix()*c;
               break;
            case(2):
               rS.rData() = (v.comp(0).data().array() * w.comp(1).data().array() - v.comp(1).data().array() * w.comp(0).data().array()).matrix()*c;
               break;
         }
      } else
      {
         switch(static_cast<int>(id))
         {
            case(0):
               rS.rData() = v.comp(1).data().array() * w.comp(2).data().array() - v.comp(2).data().array() * w.comp(1).data().array();
               break;
            case(1):
               rS.rData() = v.comp(2).data().array() * w.comp(0).data().array() - v.comp(0).data().array() * w.comp(2).data().array();
               break;
            case(2):
               rS.rData() = v.comp(0).data().array() * w.comp(1).data().array() - v.comp(1).data().array() * w.comp(0).data().array();
               break;
         }
      }
   }

   template <int COMPONENTS> void Cross::add(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         switch(static_cast<int>(id))
         {
            case(0):
               rS.rData() += (v.comp(1).data().array() * w.comp(2).data().array() - v.comp(2).data().array() * w.comp(1).data().array()).matrix()*c;
               break;
            case(1):
               rS.rData() += (v.comp(2).data().array() * w.comp(0).data().array() - v.comp(0).data().array() * w.comp(2).data().array()).matrix()*c;
               break;
            case(2):
               rS.rData() += (v.comp(0).data().array() * w.comp(1).data().array() - v.comp(1).data().array() * w.comp(0).data().array()).matrix()*c;
               break;
         }
      } else
      {
         switch(static_cast<int>(id))
         {
            case(0):
               rS.rData() += v.comp(1).data().array() * w.comp(2).data().array() - v.comp(2).data().array() * w.comp(1).data().array();
               break;
            case(1):
               rS.rData() += v.comp(2).data().array() * w.comp(0).data().array() - v.comp(0).data().array() * w.comp(2).data().array();
               break;
            case(2):
               rS.rData() += v.comp(0).data().array() * w.comp(1).data().array() - v.comp(1).data().array() * w.comp(0).data().array();
               break;
         }
      }
   }

   template <int COMPONENTS> void Cross::sub(FieldComponents::Physical::Id id, Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         switch(static_cast<int>(id))
         {
            case(0):
               rS.rData() -= (v.comp(1).data().array() * w.comp(2).data().array() - v.comp(2).data().array() * w.comp(1).data().array()).matrix()*c;
               break;
            case(1):
               rS.rData() -= (v.comp(2).data().array() * w.comp(0).data().array() - v.comp(0).data().array() * w.comp(2).data().array()).matrix()*c;
               break;
            case(2):
               rS.rData() -= (v.comp(0).data().array() * w.comp(1).data().array() - v.comp(1).data().array() * w.comp(0).data().array()).matrix()*c;
               break;
         }
      } else
      {
         switch(static_cast<int>(id))
         {
            case(0):
               rS.rData() -= v.comp(1).data().array() * w.comp(2).data().array() - v.comp(2).data().array() * w.comp(1).data().array();
               break;
            case(1):
               rS.rData() -= v.comp(2).data().array() * w.comp(0).data().array() - v.comp(0).data().array() * w.comp(2).data().array();
               break;
            case(2):
               rS.rData() -= v.comp(0).data().array() * w.comp(1).data().array() - v.comp(1).data().array() * w.comp(0).data().array();
               break;
         }
      }
   }
}
}

#endif // CROSS_HPP
