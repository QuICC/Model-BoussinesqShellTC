/** 
 * @file AssociatedLegendreFlyTransform.cpp
 * @brief Source of the implementation of the associated Legendre on-the-fly transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "PolynomialTransforms/AssociatedLegendreFlyTransform.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "Quadratures/LegendreRule.hpp"
#include "PolynomialTransforms/AssociatedLegendrePolynomial.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array AssociatedLegendreFlyTransform::generateGrid(const int size)
   {
      // Initialise grid storage
      Array grid(size);
      Array tmp(size);

      LegendreRule::computeQuadrature(grid, tmp, size);

      grid = grid.array().acos();

      return grid;
   }

   AssociatedLegendreFlyTransform::AssociatedLegendreFlyTransform()
      : mcMaxOpCols(64), mOpCols(2)
   {
   }

   AssociatedLegendreFlyTransform::~AssociatedLegendreFlyTransform()
   {
   }

   void AssociatedLegendreFlyTransform::init(AssociatedLegendreFlyTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Initialise the quadrature grid and weights
      this->initOperators();
   }

   void AssociatedLegendreFlyTransform::requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const
   {
      //
      // No possible options
      //
   }

   void AssociatedLegendreFlyTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId)
   {
      //
      // No possible options
      //
   }

   const Array& AssociatedLegendreFlyTransform::meshGrid() const
   {
      if(this->mXGrid.size() == 0 || this->mThGrid.size() == 0 || this->mWeights.size() == 0)
      {
         throw Exception("AssociatedLegendreFlyTransform has not been initialised!");
      }

      return this->mThGrid;
   }

   void AssociatedLegendreFlyTransform::initOperators()
   {
      this->mXGrid.resize(this->mspSetup->fwdSize());
      this->mThGrid.resize(this->mspSetup->fwdSize());
      this->mWeights.resize(this->mspSetup->fwdSize());

      internal::Array igrid, iweights;

      LegendreRule::computeQuadrature(this->mXGrid, this->mWeights, igrid, iweights, this->mspSetup->fwdSize());

      this->mThGrid = this->mXGrid.array().acos();

      // Normalise weights by 2*pi for spherical harmonics
      this->mWeights.array() *= 2.0*Math::PI;

      // Initialise storage for l(l+1) factor
      this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0, this->mspSetup->specSize()-1);
      this->mLl = this->mLl.array()*(this->mLl.array() + 1.0);
      this->mLl2 = this->mLl.array().pow(2);
      this->mDivLl = this->mLl.array().pow(-1);
      this->mDivLl(0) = 0.0;

      
      // Initialize recurrence storage
      int maxCols = 2;
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         maxCols =  std::max(maxCols, static_cast<int>(this->mspSetup->fast().at(iM).size()));
      }
      #ifdef GEOMHDISCC_MEMORYUSAGE_HIGH
         this->mOpCols = maxCols;
      #else
         this->mOpCols = std::min(maxCols, this->mcMaxOpCols);
      #endif //GEOMHDISCC_MEMORYUSAGE_HIGH
      this->mOp.resize(this->mspSetup->fwdSize(), this->mOpCols);
   }

   void AssociatedLegendreFlyTransform::integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, AssociatedLegendreFlyTransform::IntegratorType::Id integrator, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rSpecVal.cols() == this->mspSetup->howmany());

      // Compute first derivative integration
      if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIFF)
      {
         throw Exception("NOT YET IMPLEMENTED");

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGLLDIFF)
      { 
         throw Exception("NOT YET IMPLEMENTED");

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIVLLDIFF)
      { 
         throw Exception("NOT YET IMPLEMENTED");

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGLLDIVSIN)
      { 
         throw Exception("NOT YET IMPLEMENTED");

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIVLLDIVSIN)
      { 
         throw Exception("NOT YET IMPLEMENTED");

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIVSIN)
      { 
         throw Exception("NOT YET IMPLEMENTED");

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGLL2)
      { 
         throw Exception("NOT YET IMPLEMENTED");

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGLL)
      { 
         throw Exception("NOT YET IMPLEMENTED");

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIVLL)
      { 
         throw Exception("NOT YET IMPLEMENTED");

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTG)
      {
         this->setIntegrator(rSpecVal, physVal);

      } else
      {
         throw Exception("Requested an unknown integrator");
      }
   }

   void AssociatedLegendreFlyTransform::project(MatrixZ& rPhysVal, const MatrixZ& specVal, AssociatedLegendreFlyTransform::ProjectorType::Id projector, Arithmetics::Id arithId)
   {
      assert(arithId == Arithmetics::SET);

      // assert right sizes for input  matrix
      assert(specVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative
      if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIFF)
      {
         throw Exception("NOT YET IMPLEMENTED");

      // Compute \f$l(l+1)D\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIFFLL)
      {
         throw Exception("NOT YET IMPLEMENTED");

      // Compute \f$l(l+1)/\sin\theta\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIVSIN)
      {
         throw Exception("NOT YET IMPLEMENTED");

      // Compute \f$1/\sin\theta\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIVSINLL)
      {
         throw Exception("NOT YET IMPLEMENTED");

      // Compute \f$1/\sin\theta \partial \sin\theta\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIVSINDIFFSIN)
      {
         throw Exception("DIVSINDIFFSIN not yet implemented");

      // Compute \f$l(l+1)\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::PROJLL)
      {
         throw Exception("NOT YET IMPLEMENTED");

      // Compute simple projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::PROJ)
      {
         this->setProjector(rPhysVal, specVal);

      } else
      {
         throw Exception("Requested an unknown projector");
      }
   }

   void AssociatedLegendreFlyTransform::setIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal)
   {
      // Compute integration
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 

      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int cols = this->mspSetup->mult()(iM);
         int specRows = this->mspSetup->fast().at(iM).size();

         // Main loop
         int i0 = 0;
         for(int i = 0; i < specRows/this->mOpCols; i++)
         {
            // Build operator 
            this->buildWeightedPlm(iM, i0, std::min(this->mOpCols, specRows - i0));

            // Compute quadrature integration
            rSpecVal.block(i0, start, this->mOpCols, cols) = this->mOp.transpose()*physVal.block(0, start, physRows, cols);

            i0 = (i+1)*this->mOpCols;
         }

         // Operator is smaller than storage size
         int rem = specRows%this->mOpCols;
         if(rem > 0)
         {
            // Build operator 
            this->buildWeightedPlm(iM, i0, rem);

            rSpecVal.block(i0, start, rem, cols) = this->mOp.leftCols(rem).transpose()*physVal.block(0, start, physRows, cols);
         }

         start += cols;
      }

   }

   void AssociatedLegendreFlyTransform::setProjector(MatrixZ& rPhysVal, const MatrixZ& specVal)
   {
   }

   void AssociatedLegendreFlyTransform::buildPlm(const int idx, const int start, const int nL)
   {
      // Safety assert
      assert(this->mOp.cols() >= nL);

      // Get harmonic order
      int m = this->mspSetup->slow()(idx);

      // Set P_m^m
      if(start == 0)
      {
         Polynomial::AssociatedLegendrePolynomial::Pmm(this->mOp.col(0), m, this->mXGrid);

         // Set P_{m+1}^m
         if(nL >= 1) 
         {
            Polynomial::AssociatedLegendrePolynomial::Pmm1(this->mOp.col(1), m, this->mOp.col(0), this->mXGrid);
         }
      } else
      {
         if(this->mOpCols > 2)
         {
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(0), m, m+start, this->mOp.col(this->mOpCols-1), this->mOp.col(this->mOpCols-2), this->mXGrid);
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(1), m, m+start+1, this->mOp.col(0), this->mOp.col(this->mOpCols-1), this->mXGrid);
         } else
         {
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(0), m, m+start, this->mOp.col(1), this->mOp.col(0), this->mXGrid);
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(1), m, m+start+1, this->mOp.col(0), this->mOp.col(1), this->mXGrid);
         }
      }

      for(int j = 2; j < nL; ++j)
      {
         Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(j), m, m+start+j, this->mOp.col(j-1), this->mOp.col(j-2), this->mXGrid);
      }
   }

   void AssociatedLegendreFlyTransform::buildWeightedPlm(const int idx, const int start, const int nL)
   {
      // Safety assert
      assert(this->mOp.cols() >= nL);

      // Get harmonic order
      int m = this->mspSetup->slow()(idx);

      // Set P_m^m
      if(start == 0)
      {
         Polynomial::AssociatedLegendrePolynomial::Pmm(this->mOp.col(0), m, this->mXGrid);

         // Set P_{m+1}^m
         if(nL >= 1) 
         {
            Polynomial::AssociatedLegendrePolynomial::Pmm1(this->mOp.col(1), m, this->mOp.col(0), this->mXGrid);
         }

         // Weight polynomials
         this->mOp.col(0).array() *= this->mWeights.array();
         this->mOp.col(1).array() *= this->mWeights.array();

      } else
      {
         if(this->mOpCols > 2)
         {
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(0), m, m+start, this->mOp.col(this->mOpCols-1), this->mOp.col(this->mOpCols-2), this->mXGrid);
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(1), m, m+start+1, this->mOp.col(0), this->mOp.col(this->mOpCols-1), this->mXGrid);
         } else
         {
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(0), m, m+start, this->mOp.col(1), this->mOp.col(0), this->mXGrid);
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(1), m, m+start+1, this->mOp.col(0), this->mOp.col(1), this->mXGrid);
         }
      }

      for(int j = 2; j < nL; ++j)
      {
         Polynomial::AssociatedLegendrePolynomial::Plm(this->mOp.col(j), m, m+start+j, this->mOp.col(j-1), this->mOp.col(j-2), this->mXGrid);
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat AssociatedLegendreFlyTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage for the grid and weight
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mXGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mThGrid.size();
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*this->mWeights.size();

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
