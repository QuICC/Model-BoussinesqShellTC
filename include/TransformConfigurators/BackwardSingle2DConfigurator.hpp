/** \file BackwardSingle2DConfigurator.hpp
 *  \brief This defines the backward transform second exchange single splitting operations
 */

#ifndef BACKWARDSINGLE2DCONFIGURATOR_HPP
#define BACKWARDSINGLE2DCONFIGURATOR_HPP

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

namespace EPMPhoenix {

   /**
    * \brief This class defines the backward transform second exchange single splitting operations
    */
   class BackwardSingle2DConfigurator: public BackwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splittings::Locations::Id  SplitLocation = Splittings::Locations::SECOND;

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
         BackwardSingle2DConfigurator(){};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardSingle2DConfigurator(){};

      private:
   };

   template <typename TVariable> void BackwardSingle2DConfigurator::firstStep(PhysicalNames::Id name, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Compute physical projection if required
      if(coord.needPhysical(name))
      {
         BackwardSingle2DConfigurator::firstPhysical(rVariable, coord);
      }

      // Compute physical differential projection (gradient or curl) if required
      if(coord.needPhysicalDiff(name))
      {
         BackwardSingle2DConfigurator::firstPhysicalDiff(rVariable, coord);
      }
   }

   template <typename TVariable> void BackwardSingle2DConfigurator::secondStep(PhysicalNames::Id name, TVariable& rVariable, TransformCoordinatorType& coord)
   {
   }

   template <typename TVariable> void BackwardSingle2DConfigurator::lastStep(PhysicalNames::Id name, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Compute physical projection if required
      if(coord.needPhysical(name))
      {
         BackwardSingle2DConfigurator::lastPhysical(rVariable, coord);
      }

      // Compute physical differential projection (gradient or curl) if required
      if(coord.needPhysicalDiff(name))
      {
         BackwardSingle2DConfigurator::lastPhysicalDiff(rVariable, coord);
      }
   }

   inline void BackwardSingle2DConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void BackwardSingle2DConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      coord.communicator().converter3D().setupCommunication(packs);
   }

   inline void BackwardSingle2DConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
   }

   inline void BackwardSingle2DConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
      coord.communicator().converter3D().initiateBackwardCommunication();
   }

}

#endif // BACKWARDSINGLE2DCONFIGURATOR_HPP
