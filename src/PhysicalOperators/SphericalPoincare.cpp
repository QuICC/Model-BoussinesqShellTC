/** 
 * @file SphericalPoincare.cpp
 * @brief Source of the implementation of the spherical Poincare term
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalOperators/SphericalPoincare.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

   void SphericalPoincare::set(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& rGrid, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         rS.setZero();
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            coeff = csa*rGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.setProfile(-coeff*(phGrid + t).array().cos(), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.setProfile(coeff*std::cos(thGrid(iTh_))*(phGrid + t).array().sin(), iTh, iR);
            }
         }
      }
   }

   void SphericalPoincare::add(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& rGrid, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         //
         // Zero
         // 
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.subProfile(coeff*(phGrid + t).array().cos(), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.addProfile(coeff*std::cos(thGrid(iTh_))*(phGrid + t).array().sin(), iTh, iR);
            }
         }
      }
   }

   void SphericalPoincare::sub(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& rGrid, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         // 
         // Zero
         // 
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.addProfile(coeff*(phGrid + t).array().cos(), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.subProfile(coeff*std::cos(thGrid(iTh_))*(phGrid + t).array().sin(), iTh, iR);
            }
         }
      }
   }

}
}
