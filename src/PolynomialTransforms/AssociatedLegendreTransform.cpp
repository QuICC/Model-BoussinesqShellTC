/** 
 * @file AssociatedLegendreTransform.cpp
 * @brief Source of the implementation of the associated Legendre transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
      Array tmp(size);

      LegendreRule::computeQuadrature(grid, tmp, size);

      grid = grid.array().acos();

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

   void AssociatedLegendreTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void AssociatedLegendreTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      //
      // No possible options
      //
   }

   const Array& AssociatedLegendreTransform::meshGrid() const
   {
      if(this->mXGrid.size() == 0 || this->mThGrid.size() == 0 || this->mWeights.size() == 0)
      {
         throw Exception("AssociatedLegendreTransform has not been initialised!");
      }

      return this->mThGrid;
   }

   void AssociatedLegendreTransform::initQuadrature()
   {
      this->mXGrid.resize(this->mspSetup->fwdSize());
      this->mThGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      LegendreRule::computeQuadrature(this->mXGrid, this->mWeights, this->mspSetup->fwdSize());

      this->mThGrid = this->mXGrid.array().acos();

      // Normalise weights by 2*pi
      this->mWeights.array() *= 2*Math::PI;
   }

   void AssociatedLegendreTransform::initProjector()
   {
      // Reserve storage for the projectors
      this->mProj.reserve(this->mspSetup->slow().size());

      // Loop over harmonic orders
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int m = this->mspSetup->slow()(iM);

         // Allocate memory for the projector
         this->mProj.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));

         // Loop over harmonic degrees
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            // Loop over grid points
            for(int iG = 0; iG < this->mThGrid.size(); iG++)
            {
               this->mProj.at(iM)(iG,iL) = std::tr1::sph_legendre(l,m,this->mThGrid(iG));
            }
         }
      }
   }

   void AssociatedLegendreTransform::initDerivative()
   {
      // Reserve storage for the derivative
      this->mDiff.reserve(this->mspSetup->slow().size());

      // Loop over harmonic orders
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int m = this->mspSetup->slow()(iM);

         // Allocate memory for the derivative
         this->mDiff.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));

         // Loop over harmonic degrees
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            // Loop over grid points
            for(int iG = 0; iG < this->mThGrid.size(); iG++)
            {
               if(l == 0 && m == 0)
               {
                  this->mDiff.at(iM)(iG,iL) = 0.0;
               } else if(m == 0)
               {
                  this->mDiff.at(iM)(iG,iL) = -0.5*std::sqrt(static_cast<MHDFloat>(l*(l+1)))*std::tr1::sph_legendre(l,m+1,this->mThGrid(iG));
               } else if(m == l)
               {
                  this->mDiff.at(iM)(iG,iL) = 0.5*std::sqrt(static_cast<MHDFloat>(2*l))*std::tr1::sph_legendre(l,m-1,this->mThGrid(iG));
               } else
               {
                  this->mDiff.at(iM)(iG,iL) = 0.5*(std::sqrt(static_cast<MHDFloat>((l+m)*(l-m+1)))*std::tr1::sph_legendre(l,m-1,this->mThGrid(iG)) - std::sqrt(static_cast<MHDFloat>((l-m)*(l+m+1)))*std::tr1::sph_legendre(l,m+1,this->mThGrid(iG)));
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
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mXGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mThGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mWeights.size();

      // Storage for the projector
      for(size_t i = 0; i < this->mProj.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mProj.at(i).rows()*this->mProj.at(i).cols();
      }

      // Storage for the derivative
      for(size_t i = 0; i < this->mDiff.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mDiff.at(i).rows()*this->mDiff.at(i).cols();
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
