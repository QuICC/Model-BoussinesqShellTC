/** 
 * @file IRegular3DScheme.cpp
 * @brief Source of the Chebyshev + Chebyshev scheme implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/3D/IRegular3DScheme.hpp"

// Project includes
//
#include "SpatialSchemes/Tools/RegularTools.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   const int IRegular3DScheme::DIMENSIONS = 3;

   bool IRegular3DScheme::isRegular()
   {
      return true;
   }

   IRegular3DScheme::IRegular3DScheme(const ArrayI& dim)
      : ISpatialScheme(dim.size()), mI(dim(0)), mJ(dim(1)), mK(dim(2))
   {
   }

   IRegular3DScheme::~IRegular3DScheme()
   {
   }

   void IRegular3DScheme::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, Splitting::Locations::Id flag)
   {
      // Safety assertions for default values
      assert( (id.size() == 0) || (bins.size() > 0) );
      assert( id.size() == bins.size() );
      assert( n0.size() == nN.size() );
      assert( (bins.size() == 0) || (flag != Splitting::Locations::NONE) );

      // Make sure we start with empty indexes
      fwd1D.clear();
      bwd1D.clear();
      idx2D.clear();

      // Multimap for the modes
      std::multimap<int,int> modes;

      int i0 = -1;
      int iN = -1;
      ArrayI j0, jN;
      int c0 = -1;
      int cN = -1;

      // No splitting
      if(flag == Splitting::Locations::NONE)
      {
         i0 = 0;
         iN = this->dim(transId, Dimensions::Data::DAT3D);
         j0.resize(iN);
         jN.resize(iN);
         j0.setConstant(0);
         jN.setConstant(this->dim(transId, Dimensions::Data::DAT2D));
         c0 = 0;
         cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

      // Create index list for first transform
      } else if(transId == Dimensions::Transform::TRA1D)
      {
         // Splitting is on first transform
         if(flag == Splitting::Locations::FIRST)
         {
            i0 = 0;
            iN = this->dim(transId, Dimensions::Data::DAT3D);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(0);
            jN.setConstant(this->dim(transId, Dimensions::Data::DAT2D));
            c0 = n0(0);
            cN = n0(0) + nN(0);

         // Splitting is on second transform
         } else if(flag == Splitting::Locations::SECOND)
         {
            i0 = 0;
            iN = this->dim(transId, Dimensions::Data::DAT3D);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(n0(0));
            jN.setConstant(nN(0));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         // Splitting is on both transforms
         } else if(flag == Splitting::Locations::BOTH)
         {
            i0 = n0(0);
            iN = nN(0);
            j0 = n0.tail(iN);
            jN = nN.tail(iN);
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         // Splitting is on slowest index on first transforms
         } else if(flag == Splitting::Locations::FIXED)
         {
            i0 = n0(0);
            iN = nN(0);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(0);
            jN.setConstant(this->dim(transId, Dimensions::Data::DAT2D));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*iN;
         }

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         // Splitting is on first transform
         if(flag == Splitting::Locations::FIRST || flag == Splitting::Locations::FIXED)
         {
            i0 = 0;
            iN = this->dim(transId, Dimensions::Data::DAT3D);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(n0(0));
            jN.setConstant(nN(0));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         // Splitting is on second transform
         } else if(flag == Splitting::Locations::SECOND)
         {
            i0 = n0(0);
            iN = nN(0);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(0);
            jN.setConstant(this->dim(transId, Dimensions::Data::DAT2D));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         // Splitting is on both transforms
         } else if(flag == Splitting::Locations::BOTH)
         {
            i0 = n0(0);
            iN = nN(0);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(n0(1));
            jN.setConstant(nN(1));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);
         }

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         // Splitting is on first transform
         if(flag == Splitting::Locations::FIRST || flag == Splitting::Locations::FIXED)
         {
            i0 = n0(0);
            iN = nN(0);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(0);
            jN.setConstant(this->dim(transId, Dimensions::Data::DAT2D));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         // Splitting is on second transform
         } else if(flag == Splitting::Locations::SECOND)
         {
            i0 = 0;
            iN = this->dim(transId, Dimensions::Data::DAT3D);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(0);
            jN.setConstant(this->dim(transId, Dimensions::Data::DAT2D));
            c0 = n0(0);
            cN = n0(0) + nN(0);

         // Splitting is on both transforms
         } else if(flag == Splitting::Locations::BOTH)
         {
            i0 = n0(0);
            iN = nN(0);
            j0.resize(iN);
            jN.resize(iN);
            j0 = n0.tail(iN);
            jN = nN.tail(iN);
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);
         }
      }

      // Generate map for regular indexes
      RegularTools::buildMap(modes, i0, iN, j0, jN, c0, cN);

      // Fill indexes for 2D and 3D
      RegularTools::fillIndexes2D3D(idx2D, idx3D, modes);

      // Fill indexes for 1D
      RegularTools::fillIndexes1D(fwd1D, bwd1D, idx3D, this->dim(transId, Dimensions::Data::DATF1D), this->dim(transId, Dimensions::Data::DATB1D));
   }

   int IRegular3DScheme::splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      // Splittable size for first transform splitting
      if(flag == Splitting::Locations::FIRST)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }

      // Splittable size for second transform splitting
      } else if(flag == Splitting::Locations::SECOND)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);
         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);
         }

      // Splittable size for both transforms splitting
      } else if(flag == Splitting::Locations::BOTH)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);
         }

      // Splittable size for FIXEd first transforms splitting
      } else if(flag == Splitting::Locations::FIXED)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }
      }
      
      // If none of the previous conditions were right
      throw Exception("Tried to split in a unknown dimension for 3D regular case");

      return -1;
   }
}
}
