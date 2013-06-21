/** \file AssociatedLegendreTransform.cpp
 *  \brief Source of the implementation of the associated Legendre transform
 */

// System includes
//
#include <cassert>
#include <tr1/cmath>

// External includes
//

// Class include
//
#include "PolynomialTransforms/AssociatedLegendreTransform.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "Quadratures/LegendreRule.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array AssociatedLegendreTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);

      // Create Chebyshev grid
      for(int k = 0; k < size; k++)
      {
         grid(k) = std::cos((MathConstants::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(size));
      }

      return grid;
   }

   AssociatedLegendreTransform::AssociatedLegendreTransform()
   {
   }

   AssociatedLegendreTransform::~AssociatedLegendreTransform()
   {
   }

   void AssociatedLegendreTransform::init(AssociatedLegendreTransform::SharedSetupType spSetup)
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

   void AssociatedLegendreTransform::requiredOptions(std::set<NonDimensional::Id>& list) const
   {
      //
      // No possible options
      //
   }

   void AssociatedLegendreTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options)
   {
      //
      // No possible options
      //
   }

   const Array& AssociatedLegendreTransform::meshGrid() const
   {
      if(this->mGrid.size() == 0 || this->mWeights.size() == 0)
      {
         throw Exception("AssociatedLegendreTransform has not been initialised!");
      }

      return this->mGrid;
   }

   void AssociatedLegendreTransform::initQuadrature()
   {
      this->mGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      LegendreRule::computeQuadrature(this->mGrid, this->mWeights, this->mspSetup->fwdSize());
   }

   void AssociatedLegendreTransform::initProjector()
   {
      // Reserve storage for the projectors
      this->mProjector.reserve(this->mspSetup->slow().size());

      // Loop over harmonic orders
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int m = this->mspSetup->slow()(iM);

         // Allocate memory for the projector
         this->mProjector.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iM).size()));

         // Loop over harmonic degrees
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            // Loop over grid points
            for(int iG = 0; iG < this->mGrid.size(); iG++)
            {
               this->mProjector.at(iM)(iG,iL) = std::tr1::sph_legendre(l,m,this->mGrid(iG));
            }
         }
      }
   }

   void AssociatedLegendreTransform::initDerivative()
   {
      // Reserve storage for the derivative
      this->mDerivative.reserve(this->mspSetup->slow().size());

      // Loop over harmonic orders
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int m = this->mspSetup->slow()(iM);

         // Allocate memory for the derivative
         this->mDerivative.push_back(Matrix(this->mGrid.size(), this->mspSetup->fast().at(iM).size()));

         // Loop over harmonic degrees
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            // Loop over grid points
            for(int iG = 0; iG < this->mGrid.size(); iG++)
            {
               if(l == 0 && m == 0)
               {
                  this->mDerivative.at(iM)(iG,iL) = 0.0;
               } else if(m == 0)
               {
                  this->mDerivative.at(iM)(iG,iL) = -0.5*std::sqrt(static_cast<MHDFloat>(l*(l+1)))*std::tr1::sph_legendre(l,m+1,this->mGrid(iG));
               } else if(m == l)
               {
                  this->mDerivative.at(iM)(iG,iL) = 0.5*std::sqrt(static_cast<MHDFloat>(2*l))*std::tr1::sph_legendre(l,m-1,this->mGrid(iG));
               } else
               {
                  this->mDerivative.at(iM)(iG,iL) = 0.5*(std::sqrt(static_cast<MHDFloat>((l+m)*(l-m+1)))*std::tr1::sph_legendre(l,m-1,this->mGrid(iG)) - std::sqrt(static_cast<MHDFloat>((l-m)*(l+m+1)))*std::tr1::sph_legendre(l,m+1,this->mGrid(iG)));
               }
            }
         }
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat AssociatedLegendreTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage for the grid and weight
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mWeights.size();

      // Storage for the projector
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mProjector.size();

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
