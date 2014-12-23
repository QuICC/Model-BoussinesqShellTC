/** 
 * @file AssociatedLegendreTransform.hpp
 * @brief Implementation of the associated Legendre transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ASSOCIATEDLEGENDRETRANSFORM_HPP
#define ASSOCIATEDLEGENDRETRANSFORM_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//

// System includes
//
#include <set>
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"
#include "Enums/Arithmetics.hpp"
#include "Enums/NonDimensional.hpp"
#include "PolynomialTransforms/PolySetup.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Simple struct holding details about ChebyshevFFT transform
    */
   struct AssociatedLegendreIds {

      /**
       * @brief Simple struct holding the projector IDs
       */
      struct Projectors
      {
         /// Enum of projector IDs
         enum Id {PROJ, DIFF, DIVSIN, DIVSINDIFFSIN};
      };

      /**
       * @brief Simple struct holding the integrator IDs
       */
      struct Integrators
      {
         /// Enum of integrator IDs
         enum Id {INTG};
      };

   };

   /**
    * @brief Implementation of the FFTW transform for a Chebyshev expansion
    */ 
   class AssociatedLegendreTransform
   {
      public:
         /// Typedef for the configuration class
         typedef PolySetup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedPolySetup SharedSetupType;

         /// Typedef for the Projector type
         typedef AssociatedLegendreIds::Projectors ProjectorType;

         /// Typedef for the Integrator type
         typedef AssociatedLegendreIds::Integrators IntegratorType;

         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size); 

         /**
          * @brief Constructor
          */
         AssociatedLegendreTransform();

         /**
          * @brief Destructor
          */
         ~AssociatedLegendreTransform();

         /**
          * @brief Initialise the polynomial transform (matrices, weights, grid, etc)
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief set list of required options
          */
         void requiredOptions(std::set<NonDimensional::Id>& list, const Dimensions::Transform::Id dimId) const;

         /**
          * @brief Set the required options
          */
         void setOptions(const std::map<NonDimensional::Id, MHDFloat>& options, const Dimensions::Transform::Id dimId);

         /**
          * @brief Get the physical grid
          */
         const Array& meshGrid() const; 

         /**
          * @brief Compute quadrature integration
          *
          * @param rSpecVal   Output spectral coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, IntegratorType::Id integrator);

         /**
          * @brief Compute polynomial projection
          *
          * @param rPhysVal   Output physical values
          * @param specVal    Input spectral coefficients
          * @param projector  Projector to use
          *
          * @tparam TOperation   Arithmetic operation to perform
          */
         template <Arithmetics::Id TOperation> void project(MatrixZ& rPhysVal, const MatrixZ& specVal, ProjectorType::Id projector);

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief Initialise the operators
          */
         void initOperators();

         /**
          * @brief Storage for the quadrature points x = [-1, 1]
          */
         Array mXGrid;

         /**
          * @brief Storage for the quadrature points th = [0, pi]
          */
         Array mThGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         Array mWeights;

         /**
          * @brief Polynomial setup object providing the sizes
          */
         SharedPolySetup    mspSetup;

         /**
          * @brief Projector matrix
          */
         std::vector<Matrix>  mProj;

         /**
          * @brief Derivative matrix
          */
         std::vector<Matrix>  mDiff;

         /**
          * @brief \f$1/\sin\theta\f$ projector matrix
          */
         std::vector<Matrix>  mDivSin;
   };

   template <Arithmetics::Id TOperation> void AssociatedLegendreTransform::integrate(MatrixZ& rSpecVal, const MatrixZ& physVal, AssociatedLegendreTransform::IntegratorType::Id integrator)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // assert right sizes for input matrix
      assert(physVal.rows() == this->mspSetup->fwdSize());
      assert(physVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rSpecVal.cols() == this->mspSetup->howmany());

      // Compute integration
      int start = 0;
      int physRows = this->mspSetup->fwdSize(); 
      for(size_t i = 0; i < this->mProj.size(); i++)
      {
         int cols = this->mspSetup->mult()(i);
         int specRows = this->mProj.at(i).cols();
         rSpecVal.block(0, start, specRows, cols) = this->mProj.at(i).transpose()*this->mWeights.asDiagonal()*physVal.block(0,start, physRows, cols);
         start += cols;
      }
   }

   template <Arithmetics::Id TOperation> void AssociatedLegendreTransform::project(MatrixZ& rPhysVal, const MatrixZ& specVal, AssociatedLegendreTransform::ProjectorType::Id projector)
   {
      // Add static assert to make sure only SET operation is used
      Debug::StaticAssert< (TOperation == Arithmetics::SET) >();

      // assert right sizes for input  matrix
      assert(specVal.cols() == this->mspSetup->howmany());

      // assert right sizes for output matrix
      assert(rPhysVal.rows() == this->mspSetup->fwdSize());
      assert(rPhysVal.cols() == this->mspSetup->howmany());

      // Compute first derivative
      if(projector == AssociatedLegendreTransform::ProjectorType::DIFF)
      {
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDiff.at(i).cols();
            rPhysVal.block(0, start, physRows, cols) = this->mDiff.at(i)*specVal.block(0,start, specRows, cols);
            start += cols;
         }

      // Compute \f$1/\sin\theta\f$ projection
      } else if(projector == AssociatedLegendreTransform::ProjectorType::DIVSIN)
      {
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mDiff.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mDiff.at(i).cols();
            rPhysVal.block(0, start, physRows, cols) = this->mDivSin.at(i)*specVal.block(0,start, specRows, cols);
            start += cols;
         }

      // Compute simple projection
      } else
      {
         int start = 0;
         int physRows = this->mspSetup->fwdSize(); 
         for(size_t i = 0; i < this->mProj.size(); i++)
         {
            int cols = this->mspSetup->mult()(i);
            int specRows = this->mProj.at(i).cols();
            rPhysVal.block(0, start, physRows, cols) = this->mProj.at(i)*specVal.block(0,start, specRows, cols);
            start += cols;
         }
      }
   }

}
}

#endif // ASSOCIATEDLEGENDRETRANSFORM_HPP
