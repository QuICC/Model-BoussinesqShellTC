/** 
 * @file ScalarSelector.hpp
 * @brief Small template to select the right scalar field types
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

      // Configure code to use TTT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TTT
         template<> struct ScalarSelector<Dimensions::Transform::TRA1D>
         {
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> FwdType;

            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> BwdType;
         };

         template<> struct ScalarSelector<Dimensions::Transform::TRA2D>
         {
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> FwdType;

            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> BwdType;
         };

         template<> struct ScalarSelector<Dimensions::Transform::TRA3D>
         {
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> FwdType;
   
            typedef  FlatScalarField<MHDFloat, Dimensions::THREED> BwdType;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_TTT

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

      // Configure code to use FFF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_FFF
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
      #endif //GEOMHDISCC_SPATIALSCHEME_FFF

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

      // Configure code to use AFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_AFT
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
      #endif //GEOMHDISCC_SPATIALSCHEME_AFT

      // Configure code to use BLF scheme
      #if defined GEOMHDISCC_SPATIALSCHEME_BLFL || defined GEOMHDISCC_SPATIALSCHEME_BLFM
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
      #endif //GEOMHDISCC_SPATIALSCHEME_BLF

      // Configure code to use SLF scheme
      #if defined GEOMHDISCC_SPATIALSCHEME_SLFL || defined GEOMHDISCC_SPATIALSCHEME_SLFM
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
      #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL || defined GEOMHDISCC_SPATIALSCHEME_SLFM

      // Configure code to use WFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_WFT
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
      #endif //GEOMHDISCC_SPATIALSCHEME_WFT

      // Configure code to use WLF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_WLF
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
      #endif //GEOMHDISCC_SPATIALSCHEME_WLF

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
