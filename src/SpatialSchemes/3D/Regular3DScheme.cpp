/** \file Regular3DScheme.cpp
 *  \brief Source of the Chebyshev + Chebyshev scheme implementation
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/3D/Regular3DScheme.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   const int Regular3DScheme::DIMENSIONS = 3;

   Regular3DScheme::Regular3DScheme(const ArrayI& dim)
      : SpatialScheme(dim.size()), mI(dim(0)), mJ(dim(1)), mK(dim(2))
   {
   }

   Regular3DScheme::~Regular3DScheme()
   {
   }

   void Regular3DScheme::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, Splitting::Locations::Id flag)
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
         iN = this->dim3D(transId);
         j0.resize(iN);
         jN.resize(iN);
         j0.setConstant(0);
         jN.setConstant(this->dim2D(transId));
         c0 = 0;
         cN = this->dim2D(transId)*this->dim3D(transId);

      // Create index list for first transform
      } else if(transId == Dimensions::Transform::TRA1D)
      {
         // Splitting is on first transform
         if(flag == Splitting::Locations::FIRST)
         {
            i0 = 0;
            iN = this->dim3D(transId);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(0);
            jN.setConstant(this->dim2D(transId));
            c0 = n0(0);
            cN = n0(0) + nN(0);

         // Splitting is on second transform
         } else if(flag == Splitting::Locations::SECOND)
         {
            i0 = 0;
            iN = this->dim3D(transId);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(n0(0));
            jN.setConstant(nN(0));
            c0 = 0;
            cN = this->dim2D(transId)*this->dim3D(transId);

         // Splitting is on both transforms
         } else if(flag == Splitting::Locations::BOTH)
         {
            i0 = n0(0);
            iN = nN(0);
            j0 = n0.tail(iN);
            jN = nN.tail(iN);
            c0 = 0;
            cN = this->dim2D(transId)*this->dim3D(transId);
         }

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         // Splitting is on first transform
         if(flag == Splitting::Locations::FIRST)
         {
            i0 = 0;
            iN = this->dim3D(transId);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(n0(0));
            jN.setConstant(nN(0));
            c0 = 0;
            cN = this->dim2D(transId)*this->dim3D(transId);

         // Splitting is on second transform
         } else if(flag == Splitting::Locations::SECOND)
         {
            i0 = n0(0);
            iN = nN(0);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(0);
            jN.setConstant(this->dim2D(transId));
            c0 = 0;
            cN = this->dim2D(transId)*this->dim3D(transId);

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
            cN = this->dim2D(transId)*this->dim3D(transId);
         }

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         // Splitting is on first transform
         if(flag == Splitting::Locations::FIRST)
         {
            i0 = n0(0);
            iN = nN(0);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(0);
            jN.setConstant(this->dim2D(transId));
            c0 = 0;
            cN = this->dim2D(transId)*this->dim3D(transId);

         // Splitting is on second transform
         } else if(flag == Splitting::Locations::SECOND)
         {
            i0 = 0;
            iN = this->dim3D(transId);
            j0.resize(iN);
            jN.resize(iN);
            j0.setConstant(0);
            jN.setConstant(this->dim2D(transId));
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
            cN = this->dim2D(transId)*this->dim3D(transId);
         }
      }
      // Counter
      int c = 0;

      // Loop over third dimension
      for(int i = 0; i < iN; i++)
      {
         // Loop over second dimension
         for(int j = 0; j < jN(i); j++)
         {
            // Check for first mode
            if(c >= c0)
            {
               if(c >= cN)
               {
                  break;
               } else
               {
                  modes.insert(std::make_pair(i0 + i,j0(i) + j));
               }
            }
            c++;
         }
         if(c >= cN)
         {
            break;
         }
      }

      // Multimap iterator
      std::multimap<int,int>::iterator mapIt;

      // Set to extract the 3D indexes
      std::set<int>  filter;

      // Loop over all modes
      for(mapIt = modes.begin(); mapIt != modes.end(); mapIt++)
      {
         filter.insert(mapIt->first);
      }

      // Set third dimension
      idx3D.resize(filter.size());

      // Make full list of index in third dimension
      std::set<int>::iterator setIt = filter.begin();
      for(int i = 0; i < idx3D.size(); i++)
      {
         idx3D(i) = *setIt;
         setIt++;
      }

      // Make full list of indexes for second dimension
      std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> mapRange;
      for(int i = 0; i < idx3D.size(); i++)
      {
         // Create storage for indexes
         idx2D.push_back(ArrayI(modes.count(idx3D(i))));

         // Get range
         mapRange = modes.equal_range(idx3D(i));

         // Loop over range
         int j = 0;
         for(mapIt = mapRange.first; mapIt != mapRange.second; mapIt++)
         {
            idx2D.at(i)(j) = mapIt->second;
            j++;
         }
      }

      // Make full list of indexes for first dimension
      for(int i = 0; i < this->dim3D(transId); i++)
      {
         // Create storage for indexes
         fwd1D.push_back(ArrayI(this->dimFwd(transId)));

         // Fill array with indexes
         for(int j = 0; j < fwd1D.at(i).size(); j++)
         {
            fwd1D.at(i)(j) = j;
         }

         // Create storage for indexes
         bwd1D.push_back(ArrayI(this->dimBwd(transId)));

         // Fill array with indexes
         for(int j = 0; j < bwd1D.at(i).size(); j++)
         {
            bwd1D.at(i)(j) = j;
         }
      }
   }

   int Regular3DScheme::splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      // Splittable size for first transform splitting
      if(flag == Splitting::Locations::FIRST)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim2D(transId)*this->dim3D(transId);

         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim2D(transId);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim3D(transId);
         }

      // Splittable size for second transform splitting
      } else if(flag == Splitting::Locations::SECOND)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim2D(transId);
         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim3D(transId);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim2D(transId)*this->dim3D(transId);
         }

      // Splittable size for both transforms splitting
      } else if(flag == Splitting::Locations::BOTH)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim3D(transId);

         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim3D(transId);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim2D(transId);
         }
      }
      
      // If none of the previous conditions were right
      assert(false);

      throw Exception("Tried to split in a unknown dimension for 3D regular case");

      return -1;
   }
}
