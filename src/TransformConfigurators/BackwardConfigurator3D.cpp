/** 
 * @file BackwardConfigurator3D.cpp
 * @brief Source of the implementation of the base backward configurator in 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Debug/DebuggerMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/BackwardConfigurator3D.hpp"

// Project includes
//
#include "ScalarFields/FieldTools.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   void BackwardConfigurator3D::project2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Project 2D", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd2DType* pInVar;
      if(edge.recoverInput())
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().recoverBwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveBackward<Dimensions::Transform::TRA2D>();
      }

      // Get output storage
      TransformCoordinatorType::CommunicatorType::Fwd2DType *pOutVar = 0;
      TransformCoordinatorType::CommunicatorType::Fwd2DType *pRecOutVar = 0;
      if(edge.recoverOutId() >= 0)
      {
         pOutVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().provideFwd();
         pRecOutVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().recoverFwd(edge.recoverOutId());
      } else
      {
         pOutVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().provideFwd();
      }

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2DTRA);

      // Compute projection transform for second dimension 
      coord.transform2D().project(pOutVar->rData(), pInVar->data(), edge.opId<TransformSelector<Dimensions::Transform::TRA2D>::Type::ProjectorType::Id>());

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2DTRA);

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().holdBwd(*pInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().freeBwd(*pInVar);
      }

      // Combine recovered output with new calculation
      if(pRecOutVar != 0)
      {
         Datatypes::FieldTools::combine(*pRecOutVar, *pOutVar, edge.combinedArithId());

         if(edge.combinedOutId() >= 0)
         {
            coord.communicator().storage<Dimensions::Transform::TRA2D>().holdFwd(*pRecOutVar, edge.combinedOutId());
         } else
         {
            coord.communicator().transferForward<Dimensions::Transform::TRA2D>(*pRecOutVar);
         }

      // Hold data for combination   
      } else if(edge.combinedOutId() >= 0)
      {
         if(edge.combinedArithId() == Arithmetics::SETNEG)
         {
            Datatypes::FieldTools::negative(*pOutVar);
         }

         coord.communicator().storage<Dimensions::Transform::TRA2D>().holdFwd(*pOutVar, edge.combinedOutId());
      }

      // Transfer calculation
      if(edge.arithId() != Arithmetics::NONE)
      {
         if(edge.arithId() == Arithmetics::SETNEG)
         {
            Datatypes::FieldTools::negative(*pOutVar);
         }

         coord.communicator().transferForward<Dimensions::Transform::TRA2D>(*pOutVar);

      } else if(edge.combinedOutId() < 0 || pRecOutVar != 0)
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().freeFwd(*pOutVar);
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2D);
   }

}
}
