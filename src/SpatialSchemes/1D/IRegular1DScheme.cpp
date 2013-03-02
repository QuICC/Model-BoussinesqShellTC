/** \file IRegular1DScheme.cpp
 *  \brief Source of the generic regulard 2D scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "SpatialSchemes/1D/IRegular1DScheme.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

   const int IRegular1DScheme::DIMENSIONS = 1;

   IRegular1DScheme::IRegular1DScheme(const ArrayI& dim)
      : ISpatialScheme(dim.size()), mI(dim(0))
   {
   }

   IRegular1DScheme::~IRegular1DScheme()
   {
   }

   void IRegular1DScheme::fillIndexes(Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, Splitting::Locations::Id flag)
   {
      // Assert for right transform (1D case)
      assert(transId == Dimensions::Transform::TRA1D);

      // Safety assertions for default values
      assert( bins.size() == n0.size() );
      assert( n0.size() == nN.size() );
      assert( (bins.size() == 0) || (flag != Splitting::Locations::NONE) );

      // Set unused third dimension
      idx3D.resize(0);

      // Clear second dimension
      idx2D.clear();

      // Make sure we start with empty indexes
      fwd1D.clear();
      bwd1D.clear();

      // Create single forward storage for indexes
      fwd1D.push_back(ArrayI(this->dim(transId, Dimensions::Data::DATF1D)));

      // Fill array with indexes
      for(int i = 0; i < fwd1D.at(0).size(); i++)
      {
         fwd1D.at(0)(i) = i;
      }

      // Create single backward storage for indexes
      bwd1D.push_back(ArrayI(this->dim(transId, Dimensions::Data::DATB1D)));

      // Fill array with indexes
      for(int i = 0; i < bwd1D.at(0).size(); i++)
      {
         bwd1D.at(0)(i) = i;
      }
   }

   int IRegular1DScheme::splittableTotal(Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      throw Exception("There is no splitting algorithm for 1D problems!");
   }

}
