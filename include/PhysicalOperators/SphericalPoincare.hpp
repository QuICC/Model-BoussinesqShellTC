/* 
 * @file SphericalPoincare.hpp
 * @brief Implementation of the spherical poincare term
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_SPHERICALPOINCARE_HPP
#define QUICC_SPHERICALPOINCARE_HPP

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
    * @brief Implementation of the spherical poincare term
    */
   class SphericalPoincare
   {
      public:
         /**
          * @brief Set S to Poincare term
          */
         static void set(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Add Poincare term to S
          */
         static void add(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Substract Poincare term from S
          */
         static void sub(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalPoincare();

         /**
          * @brief Empty destructor
          */
         ~SphericalPoincare();
   };
}
}

#endif // QUICC_SPHERICALPOINCARE_HPP
