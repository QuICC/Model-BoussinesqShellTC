/** \file ScalarField3D.hpp
 *  \brief Base for a 3D scalar field with flat data layout
 */

#ifndef SCALARFIELD3D_HPP
#define SCALARFIELD3D_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "ScalarFields/ScalarField2D.hpp"
#include "ScalarFields/ScalarFieldSetup.hpp"
#include "Enums/Dimensions.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Base for a 1D scalar field with flat data layout
    */
   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION = Dimensions::THREED> class ScalarField3D: public ScalarField2D<TData, TLayout, DIMENSION>
   {
      public:
         /// Typedef for the configuration class to use
         typedef ScalarFieldSetup SetupType;

         /// Static constant for the dimension of the scalar field
         static const Dimensions::Type FieldDimension = Dimensions::THREED;

         /**
          * @brief Constructor with simplified interface
          *
          * @param setup Setup object for the scalar field
          */
         ScalarField3D(SetupType setup);

         /**
          * @brief Destructor
          */
         virtual ~ScalarField3D();

         /**
          * @brief Get third dimension
          */
         int dim3D() const;
         
      protected:
         /**
          * @brief Constructor
          *
          * @param spDim1D    Size of the first dimension 
          * @param spDim2D    Size of the second dimension 
          * @param dim3D      Size of the third dimension 
          * @param hasStorage Flag to specify if memory shoud be initialised
          */
         ScalarField3D(SharedArrayI spDim1D, SharedArrayI spDim2D, const int dim3D, bool hasStorage = true);

         /**
          * @brief Size in the third dimension
          */
         int mDim3D;

      private:
   };

   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION> inline int ScalarField3D<TData, TLayout, DIMENSION>::dim3D() const
   {
      return this->mDim3D;
   }

   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION> ScalarField3D<TData, TLayout, DIMENSION>::ScalarField3D(ScalarField3D<TData, TLayout, DIMENSION>::SetupType setup)
      : ScalarField2D<TData, TLayout, DIMENSION>(setup.spDim1D(), setup.spDim2D(), false), mDim3D(setup.dim3D())
   {
      // Initialise field storage
      this->initStorage(*this->mspDim1D, *this->mspDim2D, this->mDim3D);
   }

   template <typename TData, template <typename, Dimensions::Type> class TLayout, Dimensions::Type DIMENSION> ScalarField3D<TData, TLayout, DIMENSION>::ScalarField3D(SharedArrayI spDim1D, SharedArrayI spDim2D, const int dim3D, bool hasStorage)
      : ScalarField2D<TData, TLayout, DIMENSION>(spDim1D, spDim2D, false), mDim3D(dim3D)
   {
      if(hasStorage)
      {
         // Initialise field storage
         this->initStorage(*this->mspDim1D, *this->mspDim2D, this->mDim3D);
      }
   }

}
}

#endif // SCALARFIELD3D_HPP
