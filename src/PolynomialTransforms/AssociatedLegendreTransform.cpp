/** 
 * @file AssociatedLegendreTransform.cpp
 * @brief Source of the implementation of the associated Legendre transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

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
#include "PolynomialTransforms/AssociatedLegendrePolynomial.hpp"

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
      this->initOperators();
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

   void AssociatedLegendreTransform::initOperators()
   {
      this->mXGrid.resize(this->mspSetup->fwdSize());
      this->mThGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      internal::Array igrid, iweights;

      LegendreRule::computeQuadrature(this->mXGrid, this->mWeights, igrid, iweights, this->mspSetup->fwdSize());

      this->mThGrid = this->mXGrid.array().acos();

      // Normalise weights by 2*pi for spherical harmonics
      this->mWeights.array() *= 2*Math::PI;

      // Reserve storage for the projectors, 1/sin projectors and derivative
      this->mProj.reserve(this->mspSetup->slow().size());
      this->mDivSin.reserve(this->mspSetup->slow().size());
      this->mDiff.reserve(this->mspSetup->slow().size());

      // Loop over harmonic orders
      Matrix op;
      internal::Matrix  itmp, ipoly;
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int m = this->mspSetup->slow()(iM);
         int maxL = this->mspSetup->fast().at(iM)(this->mspSetup->fast().at(iM).size()-1);

         // Allocate memory for the projector, 1/sin projector and derivatives
         this->mProj.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));
         this->mDivSin.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));
         this->mDiff.push_back(Matrix(this->mThGrid.size(), this->mspSetup->fast().at(iM).size()));

         op.resize(this->mThGrid.size(), maxL - m + 1);

         // Loop over harmonic degrees for projectors
         Polynomial::AssociatedLegendrePolynomial::Plm(op, ipoly, m, igrid);
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            this->mProj.at(iM).col(iL) = op.col(l - m);
         }

         // Loop over harmonic degrees for derivative
         Polynomial::AssociatedLegendrePolynomial::dPlm(op, itmp, m, ipoly, igrid);
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            this->mDiff.at(iM).col(iL) = op.col(l - m);
         }

         // Loop over harmonic degrees for 1/\sin\theta
         Polynomial::AssociatedLegendrePolynomial::sin_1Plm(op, itmp, m, igrid);
         for(int iL = 0; iL < this->mspSetup->fast().at(iM).size(); iL++)
         {
            int l = this->mspSetup->fast().at(iM)(iL);

            this->mDivSin.at(iM).col(iL) = op.col(l - m);
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

      // Storage for the 1/\sin\theta
      for(size_t i = 0; i < this->mDivSin.size(); i++)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mDivSin.at(i).rows()*this->mDivSin.at(i).cols();
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
