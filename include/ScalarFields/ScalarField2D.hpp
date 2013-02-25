/** \file ScalarField2D.hpp
 *  \brief Base for a 2D scalar field with flat data layout
 */

#ifndef SCALARFIELD2D_HPP
#define SCALARFIELD2D_HPP

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
#include "ScalarFields/ScalarField1D.hpp"
#include "ScalarFields/ScalarFieldSetup.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Base for a 1D scalar field with flat data layout
    */
   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION = Dimensions::TWOD> class ScalarField2D: public ScalarField1D<TData, TLayout, DIMENSION>
   {
      public:
         /// Typedef for the configuration class to use
         typedef ScalarFieldSetup SetupType;

         /// Static constant for the dimension of the scalar field
         static const Dimensions::Type FieldDimension = Dimensions::TWOD;

         /**
          * @brief Constructor with simplified interface
          *
          * @param setup Setup object for the scalar field 
          */
         ScalarField2D(SetupType setup);

         /**
          * @brief Destructor
          */
         virtual ~ScalarField2D();

         /**
          * @brief Get second dimension
          *
          * @param k Index of the third dimension 
          */
         int dim2D(const int k = 0) const;
         
      protected:
         /**
          * @brief Constructor
          *
          * @param spDim1D    Size in the first dimension
          * @param spDim2D    Size in the second dimension
          * @param hasStorage Flag to specify if memory shoud be initialised
          */
         ScalarField2D(SharedArrayI spDim1D, SharedArrayI spDim2D, bool hasStorage = true);

         /**
          * @brief Size of the second dimension
          */
         SharedArrayI mspDim2D;

      private:
   };

   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION> inline int ScalarField2D<TData, TLayout, DIMENSION>::dim2D(const int k) const
   {
      // Add assert on sizes array dimension
      assert(k < this->mspDim2D->size());

      return (*this->mspDim2D)(k);
   }

   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION> ScalarField2D<TData, TLayout, DIMENSION>::ScalarField2D(ScalarField2D<TData, TLayout, DIMENSION>::SetupType setup)
      : ScalarField1D<TData, TLayout, DIMENSION>(setup.spDim1D(), false), mspDim2D(setup.spDim2D())
   {
      // Initialise field storage
      this->initStorage(*this->mspDim1D, *this->mspDim2D, 0);
   }

   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION> ScalarField2D<TData, TLayout, DIMENSION>::ScalarField2D(SharedArrayI spDim1D, SharedArrayI spDim2D, bool hasStorage)
      : ScalarField1D<TData, TLayout, DIMENSION>(spDim1D, false), mspDim2D(spDim2D)
   {
      if(hasStorage)
      {
         // Initialise field storage
         this->initStorage(*this->mspDim1D, *this->mspDim2D, 0);
      }
   }

}
}

#endif // SCALARFIELD2D_HPP
