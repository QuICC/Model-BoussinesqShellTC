/** 
 * @file WorlandTransform.cpp
 * @brief Source of the implementation of the Worland transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "PolynomialTransforms/WorlandTransform.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array WorlandTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create Chebyshev grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = std::cos((MathConstants::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(size));

         grid(k) = std::sqrt((grid(k)+1.)/2);
      }

      return grid;
   }

   WorlandTransform::WorlandTransform()
   {
   }

   WorlandTransform::~WorlandTransform()
   {
   }

   void WorlandTransform::init(WorlandTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Initialise the quadrature grid and weights
      this->initQuadrature();

      // Initialise the projector
      this->initProjector();

      // Initialise the derivative
      this->initDerivative();
   }

   void WorlandTransform::requiredOptions(std::set<NonDimensional::Id>& list) const
   {
      //
      // No possible options
      //
   }

   void WorlandTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options)
   {
      //
      // No possible options
      //
   }

   const Array& WorlandTransform::meshGrid() const
   {
      if(this->mGrid.size() == 0 || this->mWeights.size() == 0)
      {
         throw Exception("WorlandTransform has not been initialised!");
      }

      return this->mGrid;
   }

   void WorlandTransform::initQuadrature()
   {
      this->mGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      this->mGrid = WorlandTransform::generateGrid(this->mspSetup->fwdSize());

      // Set the weights
      this->mWeights.setConstant(MathConstants::PI/this->mspSetup->fwdSize());
   }

   void WorlandTransform::initProjector()
   {
      // Reserve storage for the projectors
      this->mProjector.reserve(this->mspSetup->slow().size());

      // Loop over harmonic degrees
      for(int iL = 0; iL < this->mspSetup->slow().size(); iL++)
      {
         // Allocate memory for the projector
         this->mProjector.push_back(Matrix(this->mGrid.size(), this->mspSetup->specSize()));

         this->mProjector.at(iL).setConstant(0);
      }
   }

   void WorlandTransform::initDerivative()
   {
      // Reserve storage for the derivative
      this->mDerivative.reserve(this->mspSetup->slow().size());

      // Loop over harmonic degrees
      for(int iL = 0; iL < this->mspSetup->slow().size(); iL++)
      {
         // Allocate memory for the derivative
         this->mDerivative.push_back(Matrix(this->mGrid.size(), this->mspSetup->specSize()));

         this->mDerivative.at(iL).setConstant(0);
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat WorlandTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage for the grid and weight
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mWeights.size();

      // Storage for the projector
      for(size_t i = 0; i < this->mProjector.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mProjector.at(i).rows()*this->mProjector.at(i).cols();
      }

      // Storage for the derivative
      for(size_t i = 0; i < this->mDerivative.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mDerivative.at(i).rows()*this->mDerivative.at(i).cols();
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
