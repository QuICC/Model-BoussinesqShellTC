/* 
 * @file SphericalPrecession.hpp
 * @brief Implementation of the spherical Coriolis + precession term
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_SPHERICALPRECESSION_HPP
#define QUICC_SPHERICALPRECESSION_HPP

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
    * @brief Implementation of the spherical Coriolis + precession term
    */
   class SphericalPrecession
   {
      public:
         /**
          * @brief Set S to Coriolis + precession term
          */
         static void set(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis + precesion term to S
          */
         static void add(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

         /**
          * @brief Substract (Coriolis + precession) term from S
          */
         static void sub(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalPrecession();

         /**
          * @brief Empty destructor
          */
         ~SphericalPrecession();
   };
}
}

#endif // QUICC_SPHERICALPRECESSION_HPP
