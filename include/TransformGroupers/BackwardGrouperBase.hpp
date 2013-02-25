/** \file BackwardGrouperBase.hpp
 *  \brief This class defines some basic forward transform grouping tools
 */

#ifndef BACKWARDGROUPERBASE_HPP
#define BACKWARDGROUPERBASE_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Enums/PhysicalNames.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm
    */
   class BackwardGrouperBase
   {
      public:
         /**
          * @brief Setup the full backward transform structure
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         virtual void transform(std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectors, TransformCoordinatorType& coord) = 0;

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param varInfo Variable information
          */
         virtual ArrayI packs1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo) = 0;

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param varInfo Variable information
          */
         virtual ArrayI packs2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo) = 0;

         /**
          * @brief Location of the split in the configurator
          */
         Splittings::Locations::Id split;

      protected:
         /**
          * @brief Number of first exchange packets required for a scalar
          */
         const int mScalarPacks1D;

         /**
          * @brief Number of first exchange packets required for a scalar
          */
         const int mGradientPacks1D;

         /**
          * @brief Number of first exchange packets required for a vector field
          */
         const int mVectorPacks1D;

         /**
          * @brief Number of first exchange packets required for a vector curl
          */
         const int mCurlPacks1D;

         /**
          * @brief Number of second exchange packets required for a scalar
          */
         const int mScalarPacks2D;

         /**
          * @brief Number of second exchange packets required for a gradient
          */
         const int mGradientPacks2D;

         /**
          * @brief Number of second exchange packets required for a vector field
          */
         const int mVectorPacks2D;

         /**
          * @brief Number of second exchange packets required for a vector curl
          */
         const int mCurlPacks2D;

         /**
          * @brief Storage for named packet sizes for the first exchange
          */
         std::map<PhysicalNames::Id, int>  mNamedPacks1D;

         /**
          * @brief Storage for named packet sizes for the second exchange
          */
         std::map<PhysicalNames::Id, int>  mNamedPacks2D;

         /**
          * @brief Get the list of pack numbers for the first exchange
          *
          * @param varInfo Variable information
          */
         ArrayI listPacks1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Get the list of pack numbers for the second exchange
          *
          * @param varInfo Variable information
          */
         ArrayI listPacks2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Get and set the name pack numbers for the first exchange
          *
          * @param varInfo Variable information
          */
         ArrayI namePacks1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Get and set the named pack numbers for the second exchange
          *
          * @param varInfo Variable information
          */
         ArrayI namePacks2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Get the grouped pack number for the first exchange
          *
          * @param varInfo Variable information
          */
         ArrayI groupPacks1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Get the grouped pack number for the second exchange
          *
          * @param varInfo Variable information
          */
         ArrayI groupPacks2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Empty constructor
          */
         BackwardGrouperBase();

         /**
          * @brief Empty destructor
          */
         ~BackwardGrouperBase();

      private: 
   };

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef SharedPtrMacro<BackwardGrouperBase>   SharedBackwardGrouperBase;

}
}

#endif // BACKWARDGROUPERBASE_HPP
