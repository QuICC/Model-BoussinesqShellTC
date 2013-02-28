/** \file BackwardSingle1DConfigurator.hpp
 *  \brief This defines the backward transform first exchange single splitting operations
 *
 *  \mhdBug Needs test
 */

#ifndef BACKWARDSINGLE1DCONFIGURATOR_HPP
#define BACKWARDSINGLE1DCONFIGURATOR_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TransformConfigurators/BackwardConfigurator.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * \brief This class defines the backward transform first exchange single splitting operations
    */
   class BackwardSingle1DConfigurator: public BackwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::FIRST;

         /**
          * @brief Compute the first step in the backward transform
          *
          * @param name       Name of the field
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void firstStep(PhysicalNames::Id name, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Compute the second step in the backward transform
          *
          * @param name       Name of the field
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void secondStep(PhysicalNames::Id name, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Compute the last step in the backward transform
          *
          * @param name       Name of the field
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void lastStep(PhysicalNames::Id name, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Setup first exchange communication
          */
         static void setup1DCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Setup second exchange communication
          */
         static void setup2DCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Initiate first exchange communication
          */
         static void initiate1DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Initiate second exchange communication
          */
         static void initiate2DCommunication(TransformCoordinatorType& coord);

      protected:
         /**
          * @brief Compute the first step of the backward transform for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void firstPhysical(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Compute the first step of the backward (gradient) transform for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void firstPhysicalDiff(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Compute the first step of the backward transform for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void firstPhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Compute the first step of the backward (curl) transform for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void firstPhysicalDiff(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Compute the second step of the backward transform for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void secondPhysical(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Compute the second step of the backward (gradient) transform for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void secondPhysicalDiff(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Compute the second step of the backward transform for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void secondPhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Compute the second step of the backward (curl) transform for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void secondPhysicalDiff(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Compute the last step of the backward transform for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void lastPhysical(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Compute the last step of the backward (gradient) transform for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void lastPhysicalDiff(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Compute the last step of the backward transform for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void lastPhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Compute the last step of the backward (curl) transform for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void lastPhysicalDiff(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         BackwardSingle1DConfigurator(){};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardSingle1DConfigurator(){};

      private:
   };

   template <typename TVariable> void BackwardSingle1DConfigurator::firstStep(PhysicalNames::Id name, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Compute physical projection if required
      if(coord.needPhysical(name))
      {
         BackwardSingle1DConfigurator::firstPhysical(rVariable, coord);
      }

      // Compute physical differential projection (gradient or curl) if required
      if(coord.needPhysicalDiff(name))
      {
         BackwardSingle1DConfigurator::firstPhysicalDiff(rVariable, coord);
      }
   }

   template <typename TVariable> void BackwardSingle1DConfigurator::secondStep(PhysicalNames::Id name, TVariable& rVariable, TransformCoordinatorType& coord)
   {
   }

   template <typename TVariable> void BackwardSingle1DConfigurator::lastStep(PhysicalNames::Id name, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Compute physical projection if required
      if(coord.needPhysical(name))
      {
         BackwardSingle1DConfigurator::lastPhysical(rVariable, coord);
      }

      // Compute physical differential projection (gradient or curl) if required
      if(coord.needPhysicalDiff(name))
      {
         BackwardSingle1DConfigurator::lastPhysicalDiff(rVariable, coord);
      }
   }

   inline void BackwardSingle1DConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      coord.communicator().converter2D().setupCommunication(packs);
   }

   inline void BackwardSingle1DConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void BackwardSingle1DConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      coord.communicator().converter2D().initiateBackwardCommunication();
   }

   inline void BackwardSingle1DConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
   }

}
}

#endif // BACKWARDSINGLE1DCONFIGURATOR_HPP
