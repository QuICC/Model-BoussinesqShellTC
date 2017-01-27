/** 
 * @file SphericalPrecession.cpp
 * @brief Source of the implementation of the spherical Coriolis + precession term
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
#include "PhysicalOperators/SphericalPrecession.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

   void SphericalPrecession::set(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      MHDFloat cA = c*preC*std::sin(alpha);
      MHDFloat cB = c*(corC + perC*std::cos(alpha));
      int nR = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;

      MHDFloat theta;
      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               coeff = cA*std::cos(theta);
               rS.setProfile(coeff*(phGrid + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array(), iTh, iR);
               coeff = cB*std::sin(theta); 
               rS.subProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
               rS.addProfile(cA*(phGrid + t).array().sin()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               coeff = -cB*std::cos(theta);
               rS.setProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
               coeff = cA*std::sin(theta);
               rS.subProfile(coeff*(phGrid + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
               rS.subProfile(cA*(phGrid + t).array().sin()*v.comp(FieldComponents::Physical::R).profile(iTh,iR), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               coeff = cB*std::cos(theta);
               rS.setProfile(coeff*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR), iTh, iR);
               coeff = cA*std::sin(theta);
               rS.addProfile(coeff*(phGrid + t).array().cos()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR), iTh, iR);
               coeff = cA*std::cos(theta);
               rS.subProfile(coeff*(phGrid + t).array().cos()*v.comp(FieldComponents::Physical::R).profile(iTh,iR), iTh, iR);
               coeff = cB*std::sin(theta);
               rS.addProfile(coeff*v.comp(FieldComponents::Physical::R).profile(iTh,iR), iTh, iR);
            }
         }
      }
   }

   void SphericalPrecession::add(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      int nR = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;

      MHDFloat theta;
      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

            }
         }
      }
   }

   void SphericalPrecession::sub(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, SharedResolution spRes, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      int nR = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;

      MHDFloat theta;
      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(spRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

            }
         }
      }
   }

}
}
