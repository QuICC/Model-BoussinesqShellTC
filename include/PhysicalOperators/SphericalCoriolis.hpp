/** 
 * @file SphericalCoriolis.hpp
 * @brief Implementation of the spherical coriolis term
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHERICALCORIOLISHPP
#define SPHERICALCORIOLISHPP

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
#include "TypeSelectors/ScalarSelector.hpp"
#include "Resolutions/Resolution.hpp"

namespace GeoMHDiSCC {

namespace Physical {

   /**
    * @brief Implementation of the spherical coriolis term
    */
   class SphericalCoriolis
   {
      public:
         /**
          * @brief Set S to Coriolis term
          */
         static void set(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         static void add(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         static void sub(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalCoriolis();

         /**
          * @brief Empty destructor
          */
         ~SphericalCoriolis();
   };
}
}

#endif // SPHERICALCORIOLISHPP
