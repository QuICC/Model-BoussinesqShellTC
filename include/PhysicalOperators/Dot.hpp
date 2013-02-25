/** \file Dot.hpp
 *  \brief Implementation of a generic scalar product
 */

#ifndef DOT_HPP
#define DOT_HPP

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
    * @brief Implementation of a generic scalar product
    */
   class Dot
   {
      public:
         /**
          * @brief Set S to scalar product
          */
         template <int COMPONENTS> static void set(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add scalar product to S
          */
         template <int COMPONENTS> static void add(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract scalar product from S
          */
         template <int COMPONENTS> static void sub(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         Dot();

         /**
          * @brief Empty destructor
          */
         ~Dot();
   };

   template <int COMPONENTS> void Dot::set(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.rData() = c*(v.comp(0).data().array()*w.comp(0).data().array()).matrix();

         for(int i = 1; i < COMPONENTS; i++)
         {
            rS.rData() += c*(v.comp(i).data().array()*w.comp(i).data().array()).matrix();
         }
      } else
      {
         rS.rData() = (v.comp(0).data().array()*w.comp(0).data().array()).matrix();

         for(int i = 1; i < COMPONENTS; i++)
         {
            rS.rData() += (v.comp(i).data().array()*w.comp(i).data().array()).matrix();
         }
      }
   }

   template <int COMPONENTS> void Dot::add(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         for(int i = 0; i < COMPONENTS; i++)
         {
            rS.rData() += c*(v.comp(i).data().array()*w.comp(i).data().array()).matrix();
         }
      } else
      {
         for(int i = 0; i < COMPONENTS; i++)
         {
            rS.rData() += (v.comp(i).data().array()*w.comp(i).data().array()).matrix();
         }
      }
   }

   template <int COMPONENTS> void Dot::sub(Datatypes::PhysicalScalarType &rS, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &v, const VectorField<Datatypes::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         for(int i = 0; i < COMPONENTS; i++)
         {
            rS.rData() -= c*(v.comp(i).data().array()*w.comp(i).data().array()).matrix();
         }
      } else
      {
         for(int i = 0; i < COMPONENTS; i++)
         {
            rS.rData() -= (v.comp(i).data().array()*w.comp(i).data().array()).matrix();
         }
      }
   }
}
}

#endif // DOT_HPP
