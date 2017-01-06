/** 
 * @file IRegular2DScheme.cpp
 * @brief Source of the generic regulard 2D scheme implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "SpatialSchemes/2D/IRegular2DScheme.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "SpatialSchemes/Tools/RegularTools.hpp"

namespace QuICC {

namespace Schemes {

   const int IRegular2DScheme::DIMENSIONS = 2;

   bool IRegular2DScheme::isRegular()
   {
      return true;
   }

   IRegular2DScheme::IRegular2DScheme(const ArrayI& dim)
      : ISpatialScheme(dim.size()), mI(dim(0)), mJ(dim(1))
   {
   }

   IRegular2DScheme::~IRegular2DScheme()
   {
   }

   int IRegular2DScheme::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Splitting::Locations::Id flag)
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

      // No splitting
      if(flag == Splitting::Locations::NONE)
      {
         this->splitSerial(modes, transId);

      // Splitting is on first transform
      } else if(flag == Splitting::Locations::FIRST || flag == Splitting::Locations::COUPLED2D)
      {
         this->splitSingle1D(modes, n0, nN, transId);
      }

      // Fill indexes for 2D and 3D
      RegularTools::fillIndexes2D3D(idx2D, idx3D, modes);

      // Fill indexes for 1D
      RegularTools::fillIndexes1D(fwd1D, bwd1D, idx3D, this->dim(transId, Dimensions::Data::DATF1D), this->dim(transId, Dimensions::Data::DATB1D));

      // Set status (0 for success, 1 for failure)
      int status = 0;
      if(modes.size() == 0)
      {
         status = 1;
      }

      return status;
   }

   int IRegular2DScheme::splittableTotal(const  Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      if(flag == Splitting::Locations::FIRST || flag == Splitting::Locations::COUPLED2D)
      {
         return this->dim(transId, Dimensions::Data::DAT2D);
      }
      
      // If condition has not been mached
      throw Exception("Can't split in any other dimension than first dimension in 2D regular case");

      return -1;
   }

   void IRegular2DScheme::splitSerial(std::multimap<int,int>& modes, const Dimensions::Transform::Id transId)
   {
      if(transId == Dimensions::Transform::TRA3D)
      {
         throw Exception("Tried to work on third dimension in 2D case");
      }

      ArrayI j0 = ArrayI::Zero(1);
      ArrayI jN = ArrayI::Constant(1, this->dim(transId, Dimensions::Data::DAT2D));
      int c0 = 0;
      int cN = this->dim(transId, Dimensions::Data::DAT2D);

      // Generate map for regular indexes
      RegularTools::buildMap(modes, 0, 1, j0, jN, c0, cN);
   }

   void IRegular2DScheme::splitSingle1D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      ArrayI j0, jN;
      int c0 = -1;
      int cN = -1;

      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         j0 = ArrayI::Constant(1, n0(0));
         jN = ArrayI::Constant(1, nN(0));
         c0 = 0;
         cN = nN(0);

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         j0 = ArrayI::Constant(1, n0(0));
         jN = ArrayI::Constant(1, nN(0));
         c0 = 0;
         cN = nN(0);

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         throw Exception("Tried to work on third dimension in 2D case");
      }

      // Generate map for regular indexes
      RegularTools::buildMap(modes, 0, 1, j0, jN, c0, cN);
   }

}
}
