/** 
 * @file IRegularSHmScheme.cpp
 * @brief Source of the Regular basis + Spherical Harmonics scheme implementation with m spectral ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/3D/IRegularSHmScheme.hpp"

// Project includes
//
#include "SpatialSchemes/Tools/RegularTools.hpp"
#include "SpatialSchemes/Tools/SHTools.hpp"
#include "Resolutions/Tools/SHmIndexCounter.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   const int IRegularSHmScheme::DIMENSIONS = 3;

   bool IRegularSHmScheme::isRegular()
   {
      return false;
   }

   IRegularSHmScheme::IRegularSHmScheme(const ArrayI& dim)
      : ISpatialScheme(dim.size()), mI(dim(0)), mL(dim(1)), mM(dim(2))
   {
   }

   IRegularSHmScheme::~IRegularSHmScheme()
   {
   }

   void IRegularSHmScheme::addIndexCounter(SharedResolution spRes)
   {
      SharedSHmIndexCounter   spCounter(new SHmIndexCounter(spRes->sim(), spRes->cpu()));

      spRes->setIndexCounter(spCounter);
   }

   int IRegularSHmScheme::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Splitting::Locations::Id flag)
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
      }else if(flag == Splitting::Locations::FIRST)
      {
         this->splitSingle1D(modes, n0, nN, transId);

      // Splitting is on second transform
      } else if(flag == Splitting::Locations::SECOND)
      {
         this->splitSingle2D(modes, id, bins, n0, nN, transId);

      // Splitting is on both transforms
      } else if(flag == Splitting::Locations::BOTH)
      {
         this->splitTubular(modes, id, bins, n0, nN, transId);
      }

      // Fill indexes for 2D and 3D
      RegularTools::fillIndexes2D3D(idx2D, idx3D, modes);

      // Fill indexes for 1D
      if(transId == Dimensions::Transform::TRA1D || transId == Dimensions::Transform::TRA3D)
      {
         RegularTools::fillIndexes1D(fwd1D, bwd1D, idx3D, this->dim(transId, Dimensions::Data::DATF1D), this->dim(transId, Dimensions::Data::DATB1D));

      } else if(transId == Dimensions::Transform::TRA2D)
      {
         SHTools::fillIndexes1D(fwd1D, bwd1D, idx3D, this->dim(transId, Dimensions::Data::DATF1D), this->dim(transId, Dimensions::Data::DATB1D));
      }

      // Set status (0 for success, 1 for failure)
      int status = 0;
      if(modes.size() == 0)
      {
         status = 1;
      }

      return status;
   }

   int IRegularSHmScheme::splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      // Splittable size for first transform splitting
      if(flag == Splitting::Locations::FIRST)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);
            int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

            return SHTools::nHarmonics(nL, nM);

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
            return this->dim(transId, Dimensions::Data::DAT3D);

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
            return this->dim(transId, Dimensions::Data::DAT2D);

         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);
         }
      }
      
      // If none of the previous conditions were right
      throw Exception("Tried to split in a unknown dimension for regular spherical harmonics case");

      return -1;
   }

   void IRegularSHmScheme::splitSerial(std::multimap<int,int>& modes, const Dimensions::Transform::Id transId)
   {
      if(transId == Dimensions::Transform::TRA1D)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT2D);
         int nM = this->dim(transId, Dimensions::Data::DAT3D);

         // Get full list of harmonics mapped by harmonic order m
         SHTools::buildMMap(modes, nL, nM);
      } else
      {
         int k0 = 0;
         int kN = this->dim(transId, Dimensions::Data::DAT3D);
         ArrayI j0 = ArrayI::Zero(kN);
         ArrayI jN = ArrayI::Constant(kN,this->dim(transId, Dimensions::Data::DAT2D));
         int c0 = 0;
         int cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         RegularTools::buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
   }

   void IRegularSHmScheme::splitSingle1D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT3D);
         int nM = this->dim(transId, Dimensions::Data::DAT2D);

         // Get restricted sorted list of harmonics
         SHTools::buildMHSortedMap(modes, nL, nM, n0(0), nN(0));

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         int k0 = 0;
         int kN = this->dim(transId, Dimensions::Data::DAT3D);
         ArrayI j0 = ArrayI::Constant(kN, n0(0));
         ArrayI jN = ArrayI::Constant(kN, nN(0));
         int c0 = 0;
         int cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         RegularTools::buildMap(modes, k0, kN, j0, jN, c0, cN);

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         int k0 = n0(0);
         int kN = nN(0);
         ArrayI j0 = ArrayI::Zero(kN);
         ArrayI jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
         int c0 = 0;
         int cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         RegularTools::buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
   }

   void IRegularSHmScheme::splitSingle2D(std::multimap<int,int>& modes, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);
         int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

         // Get restricted list of harmonics
         std::vector<int> binOrders;
         SHTools::binMLLoad(binOrders, nL, nM, id(0), bins(0));

         // Fill with correct modes
         for(std::vector<int>::iterator vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int l = *vIt; l < this->dim(transId, Dimensions::Data::DAT2D); l++)
            {
               modes.insert(std::make_pair(*vIt, l));
            }
         }

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);
         int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

         // Get restricted list of harmonics
         std::vector<int> binOrders;
         SHTools::binMLLoad(binOrders, nL, nM, id(0), bins(0));

         // Fill with correct modes
         for(std::vector<int>::iterator vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int r = 0; r < this->dim(transId, Dimensions::Data::DAT2D); r++)
            {
               modes.insert(std::make_pair(*vIt, r));
            }
         }

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         int k0 = 0;
         int kN = this->dim(transId, Dimensions::Data::DAT3D);
         ArrayI j0 = ArrayI::Zero(kN);
         ArrayI jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
         int c0 = n0(0);
         int cN = n0(0) + nN(0);

         RegularTools::buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
   }

   void IRegularSHmScheme::splitTubular(std::multimap<int,int>& modes, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);
         int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

         // Get restricted list of harmonic orders
         std::vector<int> binOrders;
         SHTools::binMLLoad(binOrders, nL, nM, id(1), bins(1));

         // Fill with correct modes
         std::multimap<int,int> tmp;
         for(std::vector<int>::iterator vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int l = *vIt; l < this->dim(transId, Dimensions::Data::DAT2D); l++)
            {
               tmp.insert(std::make_pair(*vIt, l));
            }
         }

         int i = 0;
         for(std::multimap<int,int>::iterator mapIt = tmp.begin(); mapIt != tmp.end(); ++mapIt)
         {
            // add mode
            if(i % bins(0) == id(0))
            {
               modes.insert(*mapIt);
            }

            ++i;
         }

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);
         int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

         // Get restricted list of harmonics
         std::vector<int> binOrders;
         SHTools::binMLLoad(binOrders, nL, nM, id(1), bins(1));

         // Fill with correct modes
         for(std::vector<int>::iterator vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int r = 0; r < nN(1); r++)
            {
               modes.insert(std::make_pair(*vIt, n0(1) + r));
            }
         }

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         int k0 = n0(0);
         int kN = nN(0);
         ArrayI j0 = n0.tail(kN);
         ArrayI jN = nN.tail(kN);
         int c0 = 0;
         int cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         RegularTools::buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
   }

   void IRegularSHmScheme::buildMLMap(std::multimap<int,int>& harmonics, const int id, const int bins)
   {
      // Assert that more than one bin is available
      assert(bins > 1);

      // Make sure the list of harmonics is empty
      harmonics.clear();

      int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);
      int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

      // Reset the loads
      this->resetLoad(); // ISchemeCosts

      // Initialise the loads
      SHTools::initMLLoad(this->mLoadList, this->mLoad, this->mOptimalLoad, nL, nM, bins);

      // Combine load into pairs
      SHTools::combineMPairs(this->mLoadList, this->mRegularLoad, nM, bins);

      // Fill bins with pairs
      SHTools::fillMPairsBins(this->mLoadList, this->mRegularLoad, this->mLoad, bins);

      // Update Load sums
      this->updateLoad(bins); // ISchemeCosts

      // Fill bins with remaining modes
      SHTools::fillMRestBins(this->mLoadList, this->mRegularLoad, this->mLoad, this->mOptimalLoad, bins);

      // Update Load sums
      this->updateLoad(bins); // ISchemeCosts

      // Convert loads to harmonic orders
      SHTools::convertLoadToOrders(this->mRegularLoad, nL);

      // Extract correct bin
      std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> mapRange;
      mapRange = this->mRegularLoad.equal_range(id);

      // Loop over all orders
      std::multimap<int,int>::iterator mapIt;
      for(mapIt = mapRange.first; mapIt != mapRange.second; mapIt++)
      {
         for(int l = mapIt->second; l < nL; l++)
         {
            harmonics.insert(std::make_pair(mapIt->second, l));
         }
      }
   }
}
}
