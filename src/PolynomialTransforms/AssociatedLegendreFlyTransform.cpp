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

      // Initize temporary storage
      this->mTmpA.resize(this->mspSetup->fwdSize(), 2);
      this->mTmpB.resize(this->mspSetup->fwdSize(), 2);
   }

   void AssociatedLegendreFlyTransform::integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, AssociatedLegendreFlyTransform::IntegratorType::Id integrator)
   {
      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rSpecVal.cols() == this->mspSetup->howmany());

      // Compute first derivative integration
      if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIFF)
      {
         this->setIntegrator(rSpecVal, physVal, &AssociatedLegendreFlyTransform::buildDPlm);

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGLLDIFF)
      { 
         this->setMultLIntegrator(rSpecVal, physVal, this->mLl, &AssociatedLegendreFlyTransform::buildDPlm);

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIVLLDIFF)
      { 
         this->setMultLIntegrator(rSpecVal, physVal, this->mDivLl, &AssociatedLegendreFlyTransform::buildDPlm);

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGLLDIVSIN)
      { 
         this->setMultLIntegrator(rSpecVal, physVal, this->mLl, &AssociatedLegendreFlyTransform::buildSin_1Plm);

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIVLLDIVSIN)
      { 
         this->setMultLIntegrator(rSpecVal, physVal, this->mDivLl, &AssociatedLegendreFlyTransform::buildSin_1Plm);

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIVSIN)
      { 
         this->setIntegrator(rSpecVal, physVal, &AssociatedLegendreFlyTransform::buildSin_1Plm);

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGLL2)
      { 
         this->setMultLIntegrator(rSpecVal, physVal, this->mLl2, &AssociatedLegendreFlyTransform::buildPlm);

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGLL)
      { 
         this->setMultLIntegrator(rSpecVal, physVal, this->mLl, &AssociatedLegendreFlyTransform::buildPlm);

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTGDIVLL)
      { 
         this->setMultLIntegrator(rSpecVal, physVal, this->mDivLl, &AssociatedLegendreFlyTransform::buildPlm);

      } else if(integrator == AssociatedLegendreFlyTransform::IntegratorType::INTG)
      {
         this->setIntegrator(rSpecVal, physVal, &AssociatedLegendreFlyTransform::buildPlm);

      } else
      {
         throw Exception("Requested an unknown integrator");
      }
   }

   void AssociatedLegendreFlyTransform::project(MatrixZ& rPhysVal, const MatrixZ& specVal, AssociatedLegendreFlyTransform::ProjectorType::Id projector)
   {
      // assert right sizes for input  matrix
      assert(specVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative
      if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIFF)
      {
         this->setProjector(rPhysVal, specVal, &AssociatedLegendreFlyTransform::buildDPlm);

      // Compute \f$l(l+1)D\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIFFLL)
      {
         this->setMultLProjector(rPhysVal, specVal, this->mLl, &AssociatedLegendreFlyTransform::buildDPlm);

      // Compute \f$l(l+1)/\sin\theta\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIVSIN)
      {
         this->setProjector(rPhysVal, specVal, &AssociatedLegendreFlyTransform::buildSin_1Plm);

      // Compute \f$1/\sin\theta\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIVSINLL)
      {
         this->setMultLProjector(rPhysVal, specVal, this->mLl, &AssociatedLegendreFlyTransform::buildSin_1Plm);

      // Compute \f$1/\sin\theta \partial \sin\theta\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::DIVSINDIFFSIN)
      {
         throw Exception("DIVSINDIFFSIN not yet implemented");

      // Compute \f$l(l+1)\f$ projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::PROJLL)
      {
         this->setMultLProjector(rPhysVal, specVal, this->mLl, &AssociatedLegendreFlyTransform::buildPlm);

      // Compute simple projection
      } else if(projector == AssociatedLegendreFlyTransform::ProjectorType::PROJ)
      {
         this->setProjector(rPhysVal, specVal, &AssociatedLegendreFlyTransform::buildPlm);

      } else
      {
         throw Exception("Requested an unknown projector");
      }
   }

   void AssociatedLegendreFlyTransform::setIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, OperatorBuilder builder)
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
            builder(this, iM, i0, std::min(this->mOpCols, specRows - i0), true);

            // Compute quadrature integration
            rSpecVal.block(i0, start, this->mOpCols, cols) = this->mOp.transpose()*physVal.block(0, start, physRows, cols);

            i0 = (i+1)*this->mOpCols;
         }

         // Operator is smaller than storage size
         int rem = specRows%this->mOpCols;
         if(rem > 0)
         {
            // Build operator 
            builder(this, iM, i0, rem, true);

            rSpecVal.block(i0, start, rem, cols) = this->mOp.leftCols(rem).transpose()*physVal.block(0, start, physRows, cols);
         }

         start += cols;
      }
   }

   void AssociatedLegendreFlyTransform::setMultLIntegrator(MatrixZ& rSpecVal, const MatrixZ& physVal, const Array& mult, OperatorBuilder builder)
   {
      // Compute integration
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 

      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int cols = this->mspSetup->mult()(iM);
         int specRows = this->mspSetup->fast().at(iM).size();

         int l0 = this->mspSetup->slow()(iM);

         // Main loop
         int i0 = 0;
         for(int i = 0; i < specRows/this->mOpCols; i++)
         {
            // Build operator 
            builder(this, iM, i0, std::min(this->mOpCols, specRows - i0), true);

            // Compute quadrature integration
            rSpecVal.block(i0, start, this->mOpCols, cols) = mult.segment(l0+i0, this->mOpCols).asDiagonal()*(this->mOp.transpose()*physVal.block(0, start, physRows, cols));

            i0 = (i+1)*this->mOpCols;
         }

         // Operator is smaller than storage size
         int rem = specRows%this->mOpCols;
         if(rem > 0)
         {
            // Build operator 
            builder(this, iM, i0, rem, true);

            rSpecVal.block(i0, start, rem, cols) = mult.bottomRows(rem).asDiagonal()*(this->mOp.leftCols(rem).transpose()*physVal.block(0, start, physRows, cols));
         }

         start += cols;
      }
   }

   void AssociatedLegendreFlyTransform::setProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, OperatorBuilder builder)
   {
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 

      rPhysVal.setZero();
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int cols = this->mspSetup->mult()(iM);
         int specRows = this->mspSetup->fast().at(iM).size();

         // Main loop
         int i0 = 0;
         for(int i = 0; i < specRows/this->mOpCols; i++)
         {
            // Build operator 
            builder(this, iM, i0, std::min(this->mOpCols, specRows - i0), false);

            // Compute quadrature integration
            rPhysVal.block(0, start, physRows, cols) += this->mOp*specVal.block(i0, start, this->mOpCols, cols);

            i0 = (i+1)*this->mOpCols;
         }

         // Operator is smaller than storage size
         int rem = specRows%this->mOpCols;
         if(rem > 0)
         {
            // Build operator 
            builder(this, iM, i0, rem, false);

            rPhysVal.block(0, start, physRows, cols) += this->mOp.leftCols(rem)*specVal.block(i0, start, this->mOpCols, cols);
         }

         start += cols;
      }
   }

   void AssociatedLegendreFlyTransform::setMultLProjector(MatrixZ& rPhysVal, const MatrixZ& specVal, const Array& mult, OperatorBuilder builder)
   {
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 

      rPhysVal.setZero();
      for(int iM = 0; iM < this->mspSetup->slow().size(); iM++)
      {
         int cols = this->mspSetup->mult()(iM);
         int specRows = this->mspSetup->fast().at(iM).size();

         int l0 = this->mspSetup->slow()(iM);

         // Main loop
         int i0 = 0;
         for(int i = 0; i < specRows/this->mOpCols; i++)
         {
            // Build operator 
            builder(this, iM, i0, std::min(this->mOpCols, specRows - i0), false);

            // Compute quadrature integration
            rPhysVal.block(0, start, physRows, cols) += this->mOp*(mult.segment(l0 + i0, this->mOpCols).asDiagonal()*specVal.block(i0, start, this->mOpCols, cols));

            i0 = (i+1)*this->mOpCols;
         }

         // Operator is smaller than storage size
         int rem = specRows%this->mOpCols;
         if(rem > 0)
         {
            // Build operator 
            builder(this, iM, i0, rem, false);

            rPhysVal.block(0, start, physRows, cols) += this->mOp.leftCols(rem)*(mult.segment(l0 + i0, this->mOpCols).asDiagonal()*specVal.block(i0, start, this->mOpCols, cols));
         }

         start += cols;
      }
   }

   void AssociatedLegendreFlyTransform::buildPlm(const int idx, const int start, const int nL, const bool isWeighted)
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
            Polynomial::AssociatedLegendrePolynomial::Pm1m(this->mOp.col(1), m, this->mOp.col(0), this->mXGrid);
         }

         // Weight polynomials
         if(isWeighted)
         {
            this->mOp.col(0).array() *= this->mWeights.array();

            if(nL >= 1)
            {
               this->mOp.col(1).array() *= this->mWeights.array();
            }
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

   void AssociatedLegendreFlyTransform::buildDPlm(const int idx, const int start, const int nL, const bool isWeighted)
   {
      // Safety assert
      assert(this->mOp.cols() >= nL);

      // Get harmonic order
      int m = this->mspSetup->slow()(idx);

      // Set P_m^m
      if(start == 0)
      {
         Polynomial::AssociatedLegendrePolynomial::dPmmA(this->mOp.col(0), m, this->mXGrid);

         // Set P_{m+1}^m
         if(nL >= 1) 
         {
            Polynomial::AssociatedLegendrePolynomial::Pmm(this->mTmpA.col(0), m, this->mXGrid);
            Polynomial::AssociatedLegendrePolynomial::dPm1mA(this->mOp.col(1), m, this->mTmpA.col(0), this->mOp.col(0), this->mXGrid);
         }

         if(nL >= 2)
         {
            Polynomial::AssociatedLegendrePolynomial::Pm1m(this->mTmpA.col(1), m, this->mTmpA.col(0), this->mXGrid);
         }

         // Weight polynomials
         if(isWeighted)
         {
            this->mOp.col(0).array() *= this->mWeights.array();

            if(nL >= 1)
            {
               this->mTmpA.col(0).array() *= this->mWeights.array();
               this->mOp.col(1).array() *= this->mWeights.array();
            }

            if(nL >= 2)
            {
               this->mTmpA.col(1).array() *= this->mWeights.array();
            }
         }

      } else
      {
         // Increment polynomial
         Polynomial::AssociatedLegendrePolynomial::Plm(this->mTmpA.col(0), m, m+start-1, this->mTmpA.col(1), this->mTmpA.col(0), this->mXGrid);
         Polynomial::AssociatedLegendrePolynomial::Plm(this->mTmpA.col(1), m, m+start, this->mTmpA.col(0), this->mTmpA.col(1), this->mXGrid);

         if(this->mOpCols > 2)
         {
            // Increment derivative
            Polynomial::AssociatedLegendrePolynomial::dPlmA(this->mOp.col(0), m, m+start, this->mOp.col(this->mOpCols-1), this->mOp.col(this->mOpCols-2), this->mTmpA.col(0), this->mXGrid);
            Polynomial::AssociatedLegendrePolynomial::dPlmA(this->mOp.col(1), m, m+start+1, this->mOp.col(0), this->mOp.col(this->mOpCols-1), this->mTmpA.col(1), this->mXGrid);

         } else
         {
            // Increment derivative
            Polynomial::AssociatedLegendrePolynomial::dPlmA(this->mOp.col(0), m, m+start, this->mOp.col(1), this->mOp.col(0), this->mTmpA.col(0), this->mXGrid);
            Polynomial::AssociatedLegendrePolynomial::dPlmA(this->mOp.col(1), m, m+start+1, this->mOp.col(0), this->mOp.col(1), this->mTmpA.col(1), this->mXGrid);
         }

         // Increment polynomial and swap columns
         Polynomial::AssociatedLegendrePolynomial::Plm(this->mTmpA.col(0), m, m+start+1, this->mTmpA.col(1), this->mTmpA.col(0), this->mXGrid);
         this->mTmpA.col(0).swap(this->mTmpA.col(1));
      }

      for(int j = 2; j < nL; ++j)
      {
         // Increment derivative
         Polynomial::AssociatedLegendrePolynomial::dPlmA(this->mOp.col(j), m, m+start+j, this->mOp.col(j-1), this->mOp.col(j-2), this->mTmpA.col(1), this->mXGrid);

         if(j < nL-1)
         {
            // Increment polynomial and swap columns
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mTmpA.col(0), m, m+start+j, this->mTmpA.col(1), this->mTmpA.col(0), this->mXGrid);
            this->mTmpA.col(0).swap(this->mTmpA.col(1));
         }
      }
   }

   void AssociatedLegendreFlyTransform::buildSin_1Plm(const int idx, const int start, const int nL, const bool isWeighted)
   {
      // Safety assert
      assert(this->mOp.cols() >= nL);

      // Get harmonic order
      int m = this->mspSetup->slow()(idx);

      // Set m == 0 to zero as it only shows up in combination with \parti_\phi
      if(m == 0 && start == 0)
      {
         this->mOp.setZero();
      } else
      {
         // Main loop start
         int j0;

         // Set P_m^m
         if(start == 0)
         {
            // Initialize P_{l+1}^{m+1}
            Polynomial::AssociatedLegendrePolynomial::Pmm(this->mTmpA.col(0), m+1, this->mXGrid);

            // Initialize P_{l+1}^{m-1}
            Polynomial::AssociatedLegendrePolynomial::Pmm(this->mTmpB.col(0), m-1, this->mXGrid);
            Polynomial::AssociatedLegendrePolynomial::Pm1m(this->mTmpB.col(1), m-1, this->mTmpB.col(0), this->mXGrid);
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mTmpB.col(0), m-1, m+1, this->mTmpB.col(1), this->mTmpB.col(0), this->mXGrid);
            this->mTmpB.col(0).swap(this->mTmpB.col(1));

            // Initialize \frac{1}{\sin\theta} P_l^m
            Polynomial::AssociatedLegendrePolynomial::sin_1Plm(this->mOp.col(0), m, m, this->mTmpA.col(0), this->mTmpB.col(1));

            // Set P_{m+1}^m
            if(nL >= 1) 
            {
               // Increment P_{l+1}^{m+1}
               Polynomial::AssociatedLegendrePolynomial::Pm1m(this->mTmpA.col(1), m+1, this->mTmpA.col(0), this->mXGrid);

               // Increment P_{l+1}^{m-1}
               Polynomial::AssociatedLegendrePolynomial::Plm(this->mTmpB.col(0), m-1, m+2, this->mTmpB.col(1), this->mTmpB.col(0), this->mXGrid);
               this->mTmpB.col(0).swap(this->mTmpB.col(1));

               // Increment \frac{1}{\sin\theta} P_l^m
               Polynomial::AssociatedLegendrePolynomial::sin_1Plm(this->mOp.col(1), m, m+1, this->mTmpA.col(1), this->mTmpB.col(1));
            }

            // Weight polynomials
            if(isWeighted)
            {
               this->mOp.col(0).array() *= this->mWeights.array();

               if(nL >= 1)
               {
                  this->mOp.col(1).array() *= this->mWeights.array();
               }

               if(nL >= 2)
               {
                  this->mTmpA.col(0).array() *= this->mWeights.array();
                  this->mTmpA.col(1).array() *= this->mWeights.array();
                  this->mTmpB.col(0).array() *= this->mWeights.array();
                  this->mTmpB.col(1).array() *= this->mWeights.array();
               }
            }

            // Main loop start
            j0 = 2;

         } else
         {
            // Main loop start
            j0 = 0;
         }

         for(int j = j0; j < nL; ++j)
         {
            // Increment P_{l+1}^{m+1}
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mTmpA.col(0), m+1, m + start+j+1, this->mTmpA.col(1), this->mTmpA.col(0), this->mXGrid);
            this->mTmpA.col(0).swap(this->mTmpA.col(1));

            // Increment P_{l+1}^{m-1}
            Polynomial::AssociatedLegendrePolynomial::Plm(this->mTmpB.col(0), m-1, m + start+j+1, this->mTmpB.col(1), this->mTmpB.col(0), this->mXGrid);
            this->mTmpB.col(0).swap(this->mTmpB.col(1));

            // Increment \frac{1}{\sin\theta} P_l^m
            Polynomial::AssociatedLegendrePolynomial::sin_1Plm(this->mOp.col(j), m, m + start+j, this->mTmpA.col(1), this->mTmpB.col(1));
         }
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
