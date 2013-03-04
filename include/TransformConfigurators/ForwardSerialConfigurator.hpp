/** \file ForwardSerialConfigurator.hpp
 *  \brief This class defines the forward transform serial operations
 *
 *  \mhdBug Needs test
 */

#ifndef FORWARDSERIALCONFIGURATOR_HPP
#define FORWARDSERIALCONFIGURATOR_HPP

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
#include "TransformConfigurators/ForwardConfigurator.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the forward transform serial operations
    */
   class ForwardSerialConfigurator: public ForwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::NONE;

         /**
          * @brief First step in transform, including the nonlinear interaction for a scalar field
          *
          * @param spEquation Scalar equation
          * @param coord      Transform coordinator
          */
         static void firstStep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief First step in transform, including the nonlinear interaction for a vector field
          *
          * @param spEquation Vector equation
          * @param coord      Transform coordinator
          */
         static void firstStep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Second step in transform for a scalar field
          *
          * @param spEquation Scalar equation
          * @param coord      Transform coordinator
          */
         static void secondStep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Second step in transform for a vector field
          *
          * @param spEquation Vector equation
          * @param coord      Transform coordinator
          */
         static void secondStep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Last step in transform for a scalar field
          *
          * @param spEquation Scalar equation
          * @param coord      Transform coordinator
          */
         static void lastStep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Last step in transform for a vector field
          *
          * @param spEquation Vector equation
          * @param coord      Transform coordinator
          */
         static void lastStep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief First exchange communication setup
          */
         static void setup1DCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Second Exchange communication setup
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
          * @brief Empty constructor
          */
         ForwardSerialConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardSerialConfigurator() {};

      private:
   };

   inline void ForwardSerialConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void ForwardSerialConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void ForwardSerialConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
   }

   inline void ForwardSerialConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
   }

}
}

#endif // FORWARDSERIALCONFIGURATOR_HPP
