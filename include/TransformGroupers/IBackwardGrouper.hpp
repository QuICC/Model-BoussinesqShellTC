/** 
 * @file IBackwardGrouper.hpp
 * @brief This class defines some basic forward transform grouping tools
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IBACKWARDGROUPER_HPP
#define IBACKWARDGROUPER_HPP

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
#include "Enums/FieldIds.hpp"
#include "Variables/VariableRequirement.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm
    */
   class IBackwardGrouper
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
         virtual ArrayI packs1D(const VariableRequirement& varInfo) = 0;

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param varInfo Variable information
          */
         virtual ArrayI packs2D(const VariableRequirement& varInfo) = 0;

         /**
          * @brief Location of the split in the configurator
          */
         Splitting::Locations::Id split;

      protected:
         /**
          * @brief Number of first exchange packets required for a scalar
          */
         const int mcScalarPacks1D;

         /**
          * @brief Number of first exchange packets required for a scalar
          */
         const int mcGradientPacks1D;

         /**
          * @brief Number of first exchange packets required for a vector field
          */
         const int mcVectorPacks1D;

         /**
          * @brief Number of first exchange packets required for a vector curl
          */
         const int mcCurlPacks1D;

         /**
          * @brief Number of second exchange packets required for a scalar
          */
         const int mcScalarPacks2D;

         /**
          * @brief Number of second exchange packets required for a gradient
          */
         const int mcGradientPacks2D;

         /**
          * @brief Number of second exchange packets required for a vector field
          */
         const int mcVectorPacks2D;

         /**
          * @brief Number of second exchange packets required for a vector curl
          */
         const int mcCurlPacks2D;

         /**
          * @brief Storage for named packet sizes for the first exchange
          */
         std::map<PhysicalNames::Id, int>  mNamedPacks1D;

         /**
          * @brief Storage for named packet sizes for the second exchange
          */
         std::map<PhysicalNames::Id, int>  mNamedPacks2D;

         /**
          * @brief Get and set the name pack numbers for the first exchange
          *
          * @param varInfo Variable information
          */
         ArrayI namePacks1D(const VariableRequirement& varInfo);

         /**
          * @brief Get and set the named pack numbers for the second exchange
          *
          * @param varInfo Variable information
          */
         ArrayI namePacks2D(const VariableRequirement& varInfo);

         /**
          * @brief Get the grouped pack number for the first exchange
          *
          * @param varInfo Variable information
          */
         ArrayI groupPacks1D(const VariableRequirement& varInfo);

         /**
          * @brief Get the grouped pack number for the second exchange
          *
          * @param varInfo Variable information
          */
         ArrayI groupPacks2D(const VariableRequirement& varInfo);

         /**
          * @brief Empty constructor
          */
         IBackwardGrouper();

         /**
          * @brief Empty destructor
          */
         ~IBackwardGrouper();

      private: 
   };

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef SharedPtrMacro<IBackwardGrouper>   SharedIBackwardGrouper;

}
}

#endif // IBACKWARDGROUPER_HPP
