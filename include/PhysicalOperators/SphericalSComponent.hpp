/** 
 * @file SphericalSComponent.hpp
 * @brief Implementation of the spherical S component of a field
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

#ifndef SPHERICALSCOMPONENTHPP
#define SPHERICALSCOMPONENTHPP

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

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical Z component of a field
    */
   class SphericalSComponent
   {
      public:
         /**
          * @brief Set S to Coriolis term
          */
         static void set(Datatypes::PhysicalScalarType &rS, SharedResolution spRes, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         static void add(Datatypes::PhysicalScalarType &rS, SharedResolution spRes, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         static void sub(Datatypes::PhysicalScalarType &rS, SharedResolution spRes, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalSComponent();

         /**
          * @brief Empty destructor
          */
         ~SphericalSComponent();
   };
}
}

#endif // SphericalSComponentHPP
