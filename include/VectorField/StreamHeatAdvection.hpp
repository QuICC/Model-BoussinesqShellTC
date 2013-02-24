/** \file StreamHeatAdvection.hpp
 *  \brief Implementation of a generic streamfunction advection including conducting state
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

   /**
    * \brief Implementation of a generic streamfunction advection including conducting state
    */
   class StreamHeatAdvection
   {
      public:
         /**
          * @brief Set S to streamfunction advection product
          */
         template <int COMPONENTS> static void set(Code::PhysicalScalarType &rS, const VectorField<Code::PhysicalScalarType, COMPONENTS> &dPsi, const VectorField<Code::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add streamfunction advection product to S
          */
         template <int COMPONENTS> static void add(Code::PhysicalScalarType &rS, const VectorField<Code::PhysicalScalarType, COMPONENTS> &dPsi, const VectorField<Code::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract streamfunction advection product from S
          */
         template <int COMPONENTS> static void sub(Code::PhysicalScalarType &rS, const VectorField<Code::PhysicalScalarType, COMPONENTS> &dPsi, const VectorField<Code::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         StreamHeatAdvection();

         /**
          * @Brief Empty destructor
          */
         ~StreamHeatAdvection();
   };

   template <int COMPONENTS> void StreamHeatAdvection::set(Code::PhysicalScalarType &rS, const VectorField<Code::PhysicalScalarType, COMPONENTS> &dPsi, const VectorField<Code::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.rData() = -c*(dPsi.comp(1).data().array()*w.comp(0).data().array() + dPsi.comp(1).data().array()).matrix();
                                                                              
         rS.rData() += c*(dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      } else
      {
         rS.rData() = -(dPsi.comp(1).data().array()*w.comp(0).data().array() + dPsi.comp(1).data().array()).matrix();

         rS.rData() += (dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      }
   }

   template <int COMPONENTS> void StreamHeatAdvection::add(Code::PhysicalScalarType &rS, const VectorField<Code::PhysicalScalarType, COMPONENTS> &dPsi, const VectorField<Code::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.rData() -= c*(dPsi.comp(1).data().array()*w.comp(0).data().array() + dPsi.comp(1).data().array()).matrix();

         rS.rData() += c*(dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      } else
      {
         rS.rData() -= (dPsi.comp(1).data().array()*w.comp(0).data().array() + dPsi.comp(1).data().array()).matrix();

         rS.rData() += (dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      }
   }

   template <int COMPONENTS> void StreamHeatAdvection::sub(Code::PhysicalScalarType &rS, const VectorField<Code::PhysicalScalarType, COMPONENTS> &dPsi, const VectorField<Code::PhysicalScalarType, COMPONENTS> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.rData() += c*(dPsi.comp(1).data().array()*w.comp(0).data().array() + dPsi.comp(1).data().array()).matrix();

         rS.rData() -= c*(dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      } else
      {
         rS.rData() += (dPsi.comp(1).data().array()*w.comp(0).data().array() + dPsi.comp(1).data().array()).matrix();

         rS.rData() -= (dPsi.comp(0).data().array()*w.comp(1).data().array()).matrix();
      }
   }
}

#endif // STREAMHEATADVECTION_HPP
