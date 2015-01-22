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

namespace GeoMHDiSCC {

namespace Physical {

   /**
    * @brief Implementation of a generic scalar product
    */
   class SphericalCoriolis
   {
      public:
         /**
          * @brief Set S to Coriolis term
          */
         static void set(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, const int nR, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         static void add(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, const int nR, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         static void sub(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, const int nR, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);
         
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

   inline void SphericalCoriolis::set(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, const int nR, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      if(compId == FieldComponents::Physical::R)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.setSlice(-c*(v.comp(FieldComponents::Physical::PHI).slice(iR)*sinTheta.asDiagonal()), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.setSlice(-v.comp(FieldComponents::Physical::PHI).slice(iR)*sinTheta.asDiagonal(), iR);
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.setSlice(-c*(v.comp(FieldComponents::Physical::PHI).slice(iR)*cosTheta.asDiagonal()), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.setSlice(-v.comp(FieldComponents::Physical::PHI).slice(iR)*cosTheta.asDiagonal(), iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.setSlice(c*(v.comp(FieldComponents::Physical::R).slice(iR)*sinTheta.asDiagonal()), iR);
               rS.addSlice(c*(v.comp(FieldComponents::Physical::THETA).slice(iR)*cosTheta.asDiagonal()), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.setSlice(v.comp(FieldComponents::Physical::R).slice(iR)*sinTheta.asDiagonal(), iR);
               rS.addSlice(v.comp(FieldComponents::Physical::THETA).slice(iR)*cosTheta.asDiagonal(), iR);
            }
         }
      }
   }

   inline void SphericalCoriolis::add(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, const int nR, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      if(compId == FieldComponents::Physical::R)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.subSlice(c*(v.comp(FieldComponents::Physical::PHI).slice(iR)*sinTheta.asDiagonal()), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.subSlice(v.comp(FieldComponents::Physical::PHI).slice(iR)*sinTheta.asDiagonal(), iR);
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.subSlice(c*(v.comp(FieldComponents::Physical::PHI).slice(iR)*cosTheta.asDiagonal()), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.subSlice(v.comp(FieldComponents::Physical::PHI).slice(iR)*cosTheta.asDiagonal(), iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.addSlice(c*(v.comp(FieldComponents::Physical::R).slice(iR)*sinTheta.asDiagonal()), iR);
               rS.addSlice(c*(v.comp(FieldComponents::Physical::THETA).slice(iR)*cosTheta.asDiagonal()), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.addSlice(v.comp(FieldComponents::Physical::R).slice(iR)*sinTheta.asDiagonal(), iR);
               rS.addSlice(v.comp(FieldComponents::Physical::THETA).slice(iR)*cosTheta.asDiagonal(), iR);
            }
         }
      }
   }

   inline void SphericalCoriolis::sub(Datatypes::PhysicalScalarType &rS, FieldComponents::Physical::Id compId, const int nR, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Datatypes::PhysicalScalarType, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      if(compId == FieldComponents::Physical::R)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.addSlice(c*(v.comp(FieldComponents::Physical::PHI).slice(iR)*sinTheta.asDiagonal()), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.addSlice(v.comp(FieldComponents::Physical::PHI).slice(iR)*sinTheta.asDiagonal(), iR);
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.addSlice(c*(v.comp(FieldComponents::Physical::PHI).slice(iR)*cosTheta.asDiagonal()), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.addSlice(v.comp(FieldComponents::Physical::PHI).slice(iR)*cosTheta.asDiagonal(), iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.subSlice(c*(v.comp(FieldComponents::Physical::R).slice(iR)*sinTheta.asDiagonal()), iR);
               rS.subSlice(c*(v.comp(FieldComponents::Physical::THETA).slice(iR)*cosTheta.asDiagonal()), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               rS.subSlice(v.comp(FieldComponents::Physical::R).slice(iR)*sinTheta.asDiagonal(), iR);
               rS.subSlice(v.comp(FieldComponents::Physical::THETA).slice(iR)*cosTheta.asDiagonal(), iR);
            }
         }
      }
   }
}
}

#endif // SPHERICALCORIOLISHPP
