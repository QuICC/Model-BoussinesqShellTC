/** 
 * @file ScalarFieldSetup.hpp
 * @brief Single configuration class for the different scalar fields
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SCALARFIELDSETUP_HPP
#define SCALARFIELDSETUP_HPP

// Debug includes
//
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Single configuration class for the different scalar fields
    */
   template<Dimensions::Type DIMENSION> class ScalarFieldSetup
   {
      public:
         /**
          * @brief Constructor for 3D scalar field
          */
         ScalarFieldSetup(SharedArrayI spDim1D, SharedArrayI spDim2D, const int dim3D);

         /**
          * @brief Destructor
          */
         virtual ~ScalarFieldSetup();

         /**
          * @brief Get number of rows in data storage
          */
         int dataRows() const;

         /**
          * @brief Get number of rows in data storage
          */
         int dataCols() const;

         /**
          * @brief Get column index corresponding to the given 2D and 3D indexes
          */
         int colIdx(const int j, const int k = 0) const;

         /**
          * @brief Get index of start of block for 3D index
          */
         int blockIdx(const int k) const;

         /**
          * @brief Get number of rows of block for 3D index
          */
         int blockRows(const int k) const;

         /**
          * @brief Get number of columns of block for 3D index
          */
         int blockCols(const int k) const;

         /**
          * @brief Get number of blocks
          */
         int nBlock() const;
         
      protected:

      private:
         /**
          * @brief Array of dimensions for the first data dimension
          */
         SharedArrayI   mspDim1D;

         /**
          * @brief Array of dimensions for the first data dimension
          */
         SharedArrayI   mspDim2D;

         /**
          * @brief Array of dimensions for the first data dimension
          */
         int mDim3D;
   };

   template <Dimensions::Type DIMENSION> ScalarFieldSetup<DIMENSION>::ScalarFieldSetup(SharedArrayI spDim1D, SharedArrayI spDim2D, const int dim3D)
      : mspDim1D(spDim1D), mspDim2D(spDim2D), mDim3D(dim3D)
   {
      // Safety assertions
      assert(dim3D > 0);
      assert(this->mspDim2D->size() == this->mDim3D);
      assert(this->mspDim2D->minCoeff() > 0);
      assert(this->mspDim1D->size() == this->mDim3D);
      assert(this->mspDim1D->minCoeff() > 0);
   }

   template <Dimensions::Type DIMENSION> ScalarFieldSetup<DIMENSION>::~ScalarFieldSetup()
   {
   }

   template <Dimensions::Type DIMENSION> int ScalarFieldSetup<DIMENSION>::dataRows() const
   {
      return this->mspDim1D->maxCoeff();
   }

   template <Dimensions::Type DIMENSION> int ScalarFieldSetup<DIMENSION>::dataCols() const
   {
      return this->mspDim2D->sum();
   }

   template <Dimensions::Type DIMENSION> int ScalarFieldSetup<DIMENSION>::colIdx(const int j, const int k) const
   {
      // Only makes sense in 2D and 3D
      Debug::StaticAssert< (DIMENSION == Dimensions::TWOD) || (DIMENSION == Dimensions::THREED) >();

      // Safety assert
      assert((DIMENSION != Dimensions::TWOD) || (k == 0));

      return this->mspDim2D->head(k).sum() + j;
   }

   template <Dimensions::Type DIMENSION> int ScalarFieldSetup<DIMENSION>::blockIdx(const int k) const
   {
      // Safety assert
      assert((DIMENSION != Dimensions::TWOD) || (k == 0));

      return this->mspDim2D->head(k).sum();
   }

   template <Dimensions::Type DIMENSION> int ScalarFieldSetup<DIMENSION>::blockRows(const int k) const
   {
      // Safety assert
      assert((DIMENSION != Dimensions::TWOD) || (k == 0));

      return (*this->mspDim1D)(k);
   }

   template <Dimensions::Type DIMENSION> int ScalarFieldSetup<DIMENSION>::blockCols(const int k) const
   {
      // Safety assert
      assert((DIMENSION != Dimensions::TWOD) || (k == 0));

      return (*this->mspDim2D)(k);
   }

   template <Dimensions::Type DIMENSION> int ScalarFieldSetup<DIMENSION>::nBlock() const
   {
      // Safety assert
      assert((DIMENSION != Dimensions::TWOD) || (this->mDim3D == 1));

      return this->mDim3D;
   }

}
}

#endif // SCALARFIELDSETUP_HPP
