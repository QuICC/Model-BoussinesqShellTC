/** \file Regular2DScheme.cpp
 *  \brief Source of the generic regulard 2D scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "SpatialSchemes/2D/Regular2DScheme.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

   const int Regular2DScheme::DIMENSIONS = 2;

   Regular2DScheme::Regular2DScheme(const ArrayI& dim)
      : SpatialScheme(dim.size()), mI(dim(0)), mJ(dim(1))
   {
   }

   Regular2DScheme::~Regular2DScheme()
   {
   }

   void Regular2DScheme::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, Splitting::Locations::Id flag)
   {
      // Assert for dimension
      assert(transId == Dimensions::Transform::TRA1D || (transId == Dimensions::Transform::TRA2D));

      // Safety assertions for default values
      assert( (id.size() == 0) || (bins.size() > 0) );
      assert( id.size() == bins.size() );
      assert( n0.size() == nN.size() );
      assert( (bins.size() == 0) || (flag != Splitting::Locations::NONE) );

      // Set unused third dimension
      idx3D.resize(0);

      // Make sure we start with empty indexes
      fwd1D.clear();
      bwd1D.clear();
      idx2D.clear();

      int j0 = 0;
      int jN = this->dim2D(transId);

      if(flag == Splitting::Locations::FIRST)
      {
         j0 = n0(0);
         jN = nN(0);
      } else if(flag == Splitting::Locations::SECOND)
      {
         // There is only 1 possibility in 2D
         assert(false);
      }

      // Create single array for second dimension
      idx2D.push_back(ArrayI(jN));

      // Make full list of indexes for second dimension
      for(int j = 0; j < idx2D.at(0).size(); j++)
      {
         idx2D.at(0)(j) = j0 + j;
      }

      // Make full list of indexes for first dimension
      for(int j = 0; j < idx2D.at(0).size(); j++)
      {
         // Create storage for forward indexes
         fwd1D.push_back(ArrayI(this->dimFwd(transId)));

         // Fill array with indexes
         for(int k = 0; k < this->dimFwd(transId); k++)
         {
            fwd1D.at(j)(k) = k;
         }

         // Create storage for forward indexes
         bwd1D.push_back(ArrayI(this->dimBwd(transId)));

         // Fill array with indexes
         for(int k = 0; k < this->dimBwd(transId); k++)
         {
            bwd1D.at(j)(k) = k;
         }
      }
   }

   int Regular2DScheme::splittableTotal(const  Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      if(flag == Splitting::Locations::FIRST)
      {
         return this->dim2D(transId);
      }
      
      // If condition has not been mached

      // Second splitting not possible in 2D problem
      assert(false);

      throw Exception("Can't split in any other dimension than first dimension in 2D regular case");

      return -1;
   }

}
