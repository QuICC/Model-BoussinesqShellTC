/** \file IForwardGrouper.hpp
 *  \brief This class defines some basic forward transform grouping tools
 */

#ifndef IFORWARDGROUPER_HPP
#define IFORWARDGROUPER_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformSelector.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Equations/IScalarPEquation.hpp"
#include "Equations/IVectorPEquation.hpp"
#include "Equations/IScalarDEquation.hpp"
#include "Equations/IVectorDEquation.hpp"
#include "Variables/VariableRequirement.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm
    */
   class IForwardGrouper
   {
      public:
         /**
          * @brief Setup the full forward transform structure for prognostic equations
          *
          * @param scalEqs  Vector of prognostic scalar equations
          * @param vectEqs  Vector of prognostic vector equations
          * @param coord    Transform coord
          */
         virtual void transform(std::vector<Equations::SharedIScalarPEquation>& scalEqs, std::vector<Equations::SharedIVectorPEquation>& vectEqs, TransformCoordinatorType& coord) = 0;

         /**
          * @brief Setup the full forward transform structure for diagnostic equations
          *
          * @param scalEqs  Vector of diagnostic scalar equations
          * @param vectEqs  Vector of diagnostic vector equations
          * @param coord    Transform coord
          */
         virtual void transform(std::vector<Equations::SharedIScalarDEquation>& scalEqs, std::vector<Equations::SharedIVectorDEquation>& vectEqs, TransformCoordinatorType& coord) = 0;

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
         const int mScalarPacks1D;

         /**
          * @brief Number of second exchange packets required for a scalar
          */
         const int mScalarPacks2D;

         /**
          * @brief Number of first exchange packets required for a vector field
          */
         const int mVectorPacks1D;

         /**
          * @brief Number of second exchange packets required for a vector field
          */
         const int mVectorPacks2D;

         /**
          * @brief Get the list of pack numbers for the first exchange
          *
          * @param varInfo Variable information
          */
         ArrayI listPacks1D(const VariableRequirement& varInfo);

         /**
          * @brief Get the list of pack numbers for the second exchange
          *
          * @param varInfo Variable information
          */
         ArrayI listPacks2D(const VariableRequirement& varInfo);

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
         IForwardGrouper();

         /**
          * @brief Empty destructor
          */
         ~IForwardGrouper();

      private: 
   };

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef SharedPtrMacro<IForwardGrouper>   SharedIForwardGrouper;

}
}

#endif // IFORWARDGROUPER_HPP
