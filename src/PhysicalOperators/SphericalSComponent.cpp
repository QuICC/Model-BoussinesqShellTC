/** 
 * @file SphericalSComponent.cpp
 * @brief Source of the implementation of the spherical Z component of a field
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalOperators/SphericalSComponent.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

   void SphericalSComponent::set(Datatypes::PhysicalScalarType &rS, SharedResolution spRes, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      if(c != 1.0)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.setProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
               rS.addProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
            }
         }
      } else
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.setProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               rS.addProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
            }
         }
      }
   }

   void SphericalSComponent::add(Datatypes::PhysicalScalarType &rS, SharedResolution spRes, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      if(c != 1.0)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.addProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
               rS.addProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
            }
         }
      } else
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.addProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               rS.addProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
            }
         }
      }
   }

   void SphericalSComponent::sub(Datatypes::PhysicalScalarType &rS, SharedResolution spRes, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      if(c != 1.0)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.subProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
               rS.subProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
            }
         }
      } else
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.subProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               rS.subProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
            }
         }
      }
   }

}
}
