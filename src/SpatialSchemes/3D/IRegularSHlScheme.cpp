/** 
 * @file IRegularSHlScheme.cpp
 * @brief Source of the Regular basis + Spherical Harmonics scheme implementation with l spectral ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "SpatialSchemes/3D/IRegularSHlScheme.hpp"

// Project includes
//
#include "SpatialSchemes/Tools/RegularTools.hpp"
#include "SpatialSchemes/Tools/SHTools.hpp"
#include "Resolutions/Tools/SHlIndexCounter.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   const int IRegularSHlScheme::DIMENSIONS = 3;

   bool IRegularSHlScheme::isRegular()
   {
      return false;
   }

   IRegularSHlScheme::IRegularSHlScheme(const ArrayI& dim)
      : ISpatialScheme(dim.size()), mI(dim(0)), mL(dim(1)), mM(dim(2))
   {
   }

   IRegularSHlScheme::~IRegularSHlScheme()
   {
   }

   void IRegularSHlScheme::addIndexCounter(SharedResolution spRes)
   {
      SharedSHlIndexCounter   spCounter(new SHlIndexCounter(spRes->sim(), spRes->cpu()));

      spRes->setIndexCounter(spCounter);
   }

   void IRegularSHlScheme::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Splitting::Locations::Id flag)
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
      std::multimap<int,int>::iterator mapIt;

      // Set to extract unique indexes
      std::set<int>  filter;
      std::set<int>::iterator setIt;

      // Initialise useful variables
      int k0 = -1;
      int kN = -1;
      ArrayI j0, jN;
      int c0 = -1;
      int cN = -1;
      bool isRegular = true;

      int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);
      int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // No splitting
      if(flag == Splitting::Locations::NONE)
      {
         if(transId == Dimensions::Transform::TRA1D)
         {
            // Get full list of harmonics mapped by harmonic degree l
            SHTools::buildLMap(modes, nL, nM);

            // Indexes structure is NOT regular
            isRegular = false;
         } else
         {
            k0 = 0;
            kN = this->dim(transId, Dimensions::Data::DAT3D);
            j0 = ArrayI::Zero(kN);
            jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

            // Indexes structure is regular
            isRegular = true;
         }

      // Create index list for first transform
      } else if(transId == Dimensions::Transform::TRA1D)
      {
         // Splitting is on first transform
         if(flag == Splitting::Locations::FIRST)
         {
            // Get restricted sorted list of harmonics
            SHTools::buildLHSortedMap(modes, nL, nM, n0(0), nN(0));

            // Indexes structure is not regular
            isRegular = false;

         // Splitting is on second transform
         } else if(flag == Splitting::Locations::SECOND)
         {
            // Get restricted list of harmonics
            std::multimap<int,int> tmp;
            this->buildMLMap(tmp, id(0), bins(0));

            // invert order of the map
            std::multimap<int,int>::iterator it;
            for(it = tmp.begin(); it != tmp.end(); it++)
            {
               modes.insert(std::make_pair(it->second, it->first));
            }

            // Indexes structure is not regular
            isRegular = false;

         // Splitting is on both transforms
         } else if(flag == Splitting::Locations::BOTH)
         {
            // Get restricted list of harmonics
            this->buildMLMap(modes, id(1), bins(1));

            std::multimap<int,int> tmp;

            // invert order of the map
            for(mapIt = modes.begin(); mapIt != modes.end(); mapIt++)
            {
               tmp.insert(std::make_pair(mapIt->second, mapIt->first));
            }

            // Clear modes
            modes.clear();

            int tN = 0;
            int t0 = 0;
            for(int i = 0; i < static_cast<int>(tmp.size()); i++)
            {
               if(i % bins(0) == id(0))
               {
                  tN++;
               }
               else if(i % bins(0) < id(0))
               {
                  t0++;
               }
            }

            // invert order of the map
            mapIt = tmp.begin();
            std::advance(mapIt, t0);
            for(int i = 0; i < tN; i++)
            {
               modes.insert(*mapIt);
               mapIt++;
            }

            // Indexes structure is not regular
            isRegular = false;
         }

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         // Splitting is on first transform
         if(flag == Splitting::Locations::FIRST)
         {
            k0 = 0;
            kN = this->dim(transId, Dimensions::Data::DAT3D);
            j0 = ArrayI::Constant(kN, n0(0));
            jN = ArrayI::Constant(kN, nN(0));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

            // Indexes structure is regular
            isRegular = true;

         // Splitting is on second transform
         } else if(flag == Splitting::Locations::SECOND)
         {
            // Get restricted list of harmonics
            this->buildMLMap(modes, id(0), bins(0));

            // Loop over all modes
            for(mapIt = modes.begin(); mapIt != modes.end(); mapIt++)
            {
               filter.insert(mapIt->first);
            }

            // Clear old modes
            modes.clear();

            // Fill with correct modes
            for(setIt = filter.begin(); setIt != filter.end(); setIt++)
            {
               for(int r = 0; r < this->dim(transId, Dimensions::Data::DAT2D); r++)
               {
                  modes.insert(std::make_pair(*setIt, r));
               }
            }

            // Clear filter
            filter.clear();

            // Indexes structure is not regular
            isRegular = false;

         // Splitting is on both transforms
         } else if(flag == Splitting::Locations::BOTH)
         {
            // Get restricted list of harmonics
            this->buildMLMap(modes, id(1), bins(1));

            // Loop over all modes
            for(mapIt = modes.begin(); mapIt != modes.end(); mapIt++)
            {
               filter.insert(mapIt->first);
            }

            // Clear old modes
            modes.clear();

            // Fill with correct modes
            int i = 0;
            for(setIt = filter.begin(); setIt != filter.end(); setIt++)
            {
               for(int r = 0; r < nN(1); r++)
               {
                  modes.insert(std::make_pair(*setIt, n0(1) + r));
               }
               i++;
            }

            // Clear filter
            filter.clear();

            // Indexes structure is not regular
            isRegular = false;
         }

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         // Splitting is on first transform
         if(flag == Splitting::Locations::FIRST)
         {
            k0 = n0(0);
            kN = nN(0);
            j0 = ArrayI::Zero(kN);
            jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

            // Indexes structure is regular
            isRegular = true;

         // Splitting is on second transform
         } else if(flag == Splitting::Locations::SECOND)
         {
            k0 = 0;
            kN = this->dim(transId, Dimensions::Data::DAT3D);
            j0 = ArrayI::Zero(kN);
            jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
            c0 = n0(0);
            cN = n0(0) + nN(0);

            // Indexes structure is regular
            isRegular = true;

         // Splitting is on both transforms
         } else if(flag == Splitting::Locations::BOTH)
         {
            k0 = n0(0);
            kN = nN(0);
            j0 = n0.tail(kN);
            jN = nN.tail(kN);
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

            // Indexes structure is regular
            isRegular = true;
         }
      }
 
      // Create modes list for dimensions second and third transform
      if(isRegular)
      {
         RegularTools::buildMap(modes, k0, kN, j0, jN, c0, cN);
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
   }

   int IRegularSHlScheme::splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
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
      }
      
      // If none of the previous conditions were right
      throw Exception("Tried to split in a unknown dimension for regular spherical harmonics case");

      return -1;
   }

   void IRegularSHlScheme::buildMLMap(std::multimap<int,int>& harmonics, const int id, const int bins)
   {
      // Assert that more than one bin is available
      assert(bins > 1);

      // Make sure the list of harmonics is empty
      harmonics.clear();

      int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);
      int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

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
