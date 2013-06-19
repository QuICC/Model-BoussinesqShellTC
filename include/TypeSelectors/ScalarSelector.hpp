/** \file ScalarSelector.hpp
 *  \brief Small template to select the right scalar field types
 */

#ifndef SCALARSELECTOR_HPP
#define SCALARSELECTOR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "ScalarFields/ScalarFieldSetup.hpp"
#include "ScalarFields/FlatScalarField.hpp"

namespace GeoMHDiSCC {

   /// Namespace containing the required typedefs for a simulation
   namespace Datatypes {

      template<Dimensions::Transform::Id TId> struct ScalarSelector;

      // Configure code to use TFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFT
         template<> struct ScalarSelector<Dimensions::Transform::TRA1D>
         {
            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> FwdType;

            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> BwdType;
         };

         template<> struct ScalarSelector<Dimensions::Transform::TRA2D>
         {
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> FwdType;

            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> BwdType;
         };

         template<> struct ScalarSelector<Dimensions::Transform::TRA3D>
         {
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> FwdType;
   
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> BwdType;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_TFT

      // Configure code to use TFF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFF
         template<> struct ScalarSelector<Dimensions::Transform::TRA1D>
         {
            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> FwdType;

            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> BwdType;
         };

         template<> struct ScalarSelector<Dimensions::Transform::TRA2D>
         {
            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> FwdType;

            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> BwdType;
         };

         template<> struct ScalarSelector<Dimensions::Transform::TRA3D>
         {
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> FwdType;
   
            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> BwdType;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_TFF

      // Configure code to use CFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_CFT
         template<> struct ScalarSelector<Dimensions::Transform::TRA1D>
         {
            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> FwdType;

            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> BwdType;
         };

         template<> struct ScalarSelector<Dimensions::Transform::TRA2D>
         {
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> FwdType;

            typedef  FlatScalarField<MHDComplex, Dimensions::THREED> BwdType;
         };

         template<> struct ScalarSelector<Dimensions::Transform::TRA3D>
         {
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> FwdType;
   
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> BwdType;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_CFT

      /// Typedef for the physical space scalar
      typedef ScalarSelector<Dimensions::Transform::TRA3D>::FwdType PhysicalScalarType;

      /// Typedef for the spectral space scalar
      typedef ScalarSelector<Dimensions::Transform::TRA1D>::BwdType SpectralScalarType;

      /// Typedef for a scalar field setup
      typedef ScalarFieldSetup<Dimensions::THREED> ScalarFieldSetupType;

      /// Typedef for a scalar field setup
      typedef SharedPtrMacro<ScalarFieldSetupType> SharedScalarFieldSetupType;

   }
}

#endif // SCALARSELECTOR_HPP
