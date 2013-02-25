/** \file ScalarSelector.hpp
 *  \brief Small template to select the right scalar field types
 */

#ifndef SCALARSELECTOR_HPP
#define SCALARSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "ScalarFields/FlatLayout.hpp"
#include "ScalarFields/ScalarField3D.hpp"

namespace GeoMHDiSCC {

   /// Namespace containing the required typedefs for a cartesian simulation
   namespace Datatypes {

      template<Dimensions::Transform::Id TId> struct ScalarSelector;

      template<> struct ScalarSelector<Dimensions::Transform::TRA1D>
      {
         typedef  ScalarField3D<MHDComplex, FlatLayout> FwdType;

         typedef  ScalarField3D<MHDComplex, FlatLayout> BwdType;
      };

      template<> struct ScalarSelector<Dimensions::Transform::TRA2D>
      {
         typedef  ScalarField3D<MHDComplex, FlatLayout> FwdType;

         typedef  ScalarField3D<MHDFloat, FlatLayout> BwdType;
      };

      template<> struct ScalarSelector<Dimensions::Transform::TRA3D>
      {
         typedef  ScalarField3D<MHDFloat, FlatLayout> FwdType;

         typedef  ScalarField3D<MHDFloat, FlatLayout> BwdType;
      };


      /// Typedef for the physical space scalar
      typedef ScalarSelector<Dimensions::Transform::TRA3D>::FwdType PhysicalScalarType;

      /// Typedef for the spectral space scalar
      typedef ScalarSelector<Dimensions::Transform::TRA1D>::BwdType SpectralScalarType;

   }
}

#endif // SCALARSELECTOR_HPP
