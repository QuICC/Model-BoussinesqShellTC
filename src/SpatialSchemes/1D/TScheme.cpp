/** \file TScheme.cpp
 *  \brief Source of the Chebyshev scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "SpatialSchemes/1D/TScheme.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

   const int TScheme::DIMENSIONS = 1;

   Transform::SharedFftSetup TScheme::spSetup1D(SharedResolution spRes)
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(0)->dimFwd();
      // Get spectral size of the FFT
      int specSize = spRes->sim()->dim(Dimensions::Space::SPECTRAL,0);

      // Get number of transforms
      int howmany = 0;
      for(int i = 0; i < spRes->cpu()->dim(0)->dim3D(); i++)
      {
         howmany += spRes->cpu()->dim(0)->dim2D(i);
      }

      return Transform::SharedFftSetup(new Transform::FftSetup(size, howmany, specSize, false));
   }

   TScheme::TScheme(const ArrayI& dim, const int shift)
      : SpatialScheme(dim.size()+shift), mI(dim(0))
   {
   }

   TScheme::~TScheme()
   {
   }

   void TScheme::setDimensions()
   {
      //
      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->mDimensions.at(0)(0) = this->nX();

      // Initialise backward dimension of first transform
      this->mDimensions.at(0)(1) = this->nI();
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->mDimensions.at(0)(0) = this->nX();

      // Initialise backward dimension of first transform
      this->mDimensions.at(0)(1) = this->nI();
   }

   void TScheme::setCosts(const int shift)
   {
      // Set radial costs
      this->mCosts(0 + shift) = 1.0;
   }

   void TScheme::setScalings(const int shift)
   {
      // Set first scaling
      this->mScalings(0 + shift) = 1.0;
   }

   void TScheme::setMemory(const int shift)
   {
      // Set first memory footprint
      this->mMemory(0 + shift) = 1.0;
   }

   bool TScheme::applicable() const
   {
      return true;
   }

   Array TScheme::loadWeights()
   {
      Array weights(this->DIMENSIONS);

      weights = this->mCosts.array() * this->mScalings.array();

      return weights;
   }

   double TScheme::memoryScore(SharedResolution spRes)
   {
      #ifdef GEOMHDISCC_MEMORYUSAGE_LIMITED
         return this->mMemory.prod();
      #else
         return 1.0;
      #endif //GEOMHDISCC_MEMORYUSAGE_LIMITED
   }

   void TScheme::fillIndexes(const int dim, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, Splitting::Locations::Id flag)
   {
      // Assert for dimension
      assert(dim < 1);

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
      fwd1D.push_back(ArrayI(this->dimFwd(dim)));

      // Fill array with indexes
      for(int i = 0; i < fwd1D.at(0).size(); i++)
      {
         fwd1D.at(0)(i) = i;
      }

      // Create single backward storage for indexes
      bwd1D.push_back(ArrayI(this->dimBwd(dim)));

      // Fill array with indexes
      for(int i = 0; i < bwd1D.at(0).size(); i++)
      {
         bwd1D.at(0)(i) = i;
      }
   }

   int TScheme::splittableTotal(const int dim, Splitting::Locations::Id flag)
   {
      throw Exception("There is no splitting algorithm for 1D problems!");
   }
}
