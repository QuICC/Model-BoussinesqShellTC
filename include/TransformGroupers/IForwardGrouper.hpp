/** 
 * @file IForwardGrouper.hpp
 * @brief This class defines some basic forward transform grouping tools
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
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
          * @brief Setup the full forward transform structure for equations
          *
          * @param scalEqs  Vector of scalar equations
          * @param vectEqs  Vector of vector equations
          * @param coord    Transform coord
          */
         virtual void transform(std::vector<Equations::SharedIScalarEquation>& scalEqs, std::vector<Equations::SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord) = 0;

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param varInfo Variable information
          * @param nonInfo Nonlinear requirements
          */
         virtual ArrayI packs1D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo) = 0;

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param varInfo Variable information
          * @param nonInfo Nonlinear requirements
          */
         virtual ArrayI packs2D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo) = 0;

         /**
          * @brief Location of the split in the configurator
          */
         Splitting::Locations::Id split;

      protected:
         /**
          * @brief Constant number of first exchange packets required for a scalar
          */
         const int mcScalarPacks1D;

         /**
          * @brief Constant number of second exchange packets required for a scalar
          */
         const int mcScalarPacks2D;

         /**
          * @brief Constant number of first exchange packets required for a vector field
          */
         const int mcVectorPacks1D;

         /**
          * @brief Constant number of second exchange packets required for a vector field
          */
         const int mcVectorPacks2D;

         /**
          * @brief Actual number of first exchange packets required for a scalar
          */
         int mScalarPacks1D;

         /**
          * @brief Actula number of second exchange packets required for a scalar
          */
         int mScalarPacks2D;

         /**
          * @brief Actual number of first exchange packets required for a vector field
          */
         int mVectorPacks1D;

         /**
          * @brief Actual number of second exchange packets required for a vector field
          */
         int mVectorPacks2D;

         /**
          * @brief Get the list of pack numbers for the first exchange
          *
          * @param varInfo Variable information
          * @param nonInfo Nonlinear requirements
          */
         ArrayI listPacks1D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo);

         /**
          * @brief Get the list of pack numbers for the second exchange
          *
          * @param varInfo Variable information
          * @param nonInfo Nonlinear requirements
          */
         ArrayI listPacks2D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo);

         /**
          * @brief Get the grouped pack number for the first exchange
          *
          * @param varInfo Variable information
          * @param nonInfo Nonlinear requirements
          */
         ArrayI groupPacks1D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo);

         /**
          * @brief Get the grouped pack number for the second exchange
          *
          * @param varInfo Variable information
          * @param nonInfo Nonlinear requirements
          */
         ArrayI groupPacks2D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo);

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
