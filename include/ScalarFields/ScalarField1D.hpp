/** \file ScalarField1D.hpp
 *  \brief Base for a 1D scalar field with flat data layout
 */

#ifndef SCALARFIELD1D_HPP
#define SCALARFIELD1D_HPP

// Configuration includes
//

// System includes
//
#include <assert.h>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "ScalarField/ScalarFieldSetup.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Base for a 1D scalar field with flat data layout
    */
   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION = Dimensions::ONED> class ScalarField1D: public TLayout<TData, DIMENSION>
   {
      public:
         /// Typedef for the configuration class to use
         typedef ScalarFieldSetup SetupType;

         /// Static constant for the dimension of the scalar field
         static const Dimensions::Type FieldDimension = Dimensions::ONED;

         /**
          * @brief Constructor with simplified interface
          *
          * @param setup Setup object for the scalar field 
          */
         ScalarField1D(SetupType setup);

         /**
          * @brief Destructor
          */
         virtual ~ScalarField1D();

         /**
          * @brief Get first dimension
          *
          * @param j Index of the second dimension
          * @param k Index of the third dimension
          */
         int dim1D(const int j = 0, const int k = 0) const;
         
      protected:
         /**
          * @brief Constructor
          *
          * @param spDim1D    Size of the first dimension
          * @param hasStorage Flag to specify if memory shoud be initialised
          */
         ScalarField1D(SharedArrayI spDim1D, bool hasStorage = true);

         /**
          * @brief Dimension of the field
          */
         SharedArrayI mspDim1D;

      private:
   };

   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION> inline int ScalarField1D<TData, TLayout, DIMENSION>::dim1D(const int j, const int k) const
   {
      // Add assert on sizes matrix dimensions
      assert(k < this->mspDim1D->rows());

      return (*this->mspDim1D)(k);
   }

   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION> ScalarField1D<TData, TLayout, DIMENSION>::ScalarField1D(ScalarField1D<TData, TLayout, DIMENSION>::SetupType setup)
      : TLayout<TData, DIMENSION>(), mspDim1D(setup.spDim1D())
   {
      ArrayI  noDim(1);
      noDim(0) = 0;

      // Initialise field storage
      this->initStorage(*this->mspDim1D, noDim, 0);
   }

   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION> ScalarField1D<TData, TLayout, DIMENSION>::ScalarField1D(SharedArrayI spDim1D, bool hasStorage)
      : TLayout<TData, DIMENSION>(), mspDim1D(spDim1D)
   {
      if(hasStorage)
      {
         ArrayI  noDim;

         // Initialise field storage
         this->initStorage(*this->mspDim1D, noDim, 0);
      }
   }

}
}

#endif // SCALARFIELD1D_HPP
