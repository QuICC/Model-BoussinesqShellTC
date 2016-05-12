/** 
 * @file ForwardConfigurator3D.cpp
 * @brief Source of the implementation of the base forward configurator in 3D space
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
#include "TransformConfigurators/ForwardConfigurator3D.hpp"

// Project includes
//
#include "ScalarFields/FieldTools.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   void ForwardConfigurator3D::integrate2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 2D", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD2D);

      TransformCoordinatorType::CommunicatorType::Fwd2DType* pInVar;
      if(edge.recoverInput())
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().recoverFwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveForward<Dimensions::Transform::TRA2D>();
      }

      // Get output storage
      TransformCoordinatorType::CommunicatorType::Bwd2DType *pOutVar = 0;
      TransformCoordinatorType::CommunicatorType::Bwd2DType *pRecOutVar = 0;
      if(edge.recoverOutId() >= 0)
      {
         if(edge.combinedArithId() == Arithmetics::SET || edge.combinedArithId() == Arithmetics::SETNEG)
         {
            pOutVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().recoverBwd(edge.recoverOutId());
         } else
         {
            pOutVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().provideBwd();
            pRecOutVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().recoverBwd(edge.recoverOutId());
         }
      } else
      {
         pOutVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().provideBwd();
      }

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD2DTRA);

      // Compute integration transform of second dimension
      coord.transform2D().integrate(pOutVar->rData(), pInVar->data(), edge.opId<TransformSelector<Dimensions::Transform::TRA2D>::Type::IntegratorType::Id>(), Arithmetics::SET);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD2DTRA);

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().holdFwd(*pInVar);

      // Free temporary storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().freeFwd(*pInVar);
      }

      // Combine recovered output with new calculation
      if(pRecOutVar != 0)
      {
         Datatypes::FieldTools::combine(*pRecOutVar, *pOutVar, edge.combinedArithId());

         if(edge.combinedOutId() >= 0)
         {
            coord.communicator().storage<Dimensions::Transform::TRA2D>().holdBwd(*pRecOutVar, edge.combinedOutId());
         } else
         {
            coord.communicator().transferBackward<Dimensions::Transform::TRA2D>(*pRecOutVar);
         }
      }

      // Transfer calculation
      if(edge.arithId() != Arithmetics::NONE)
      {
         if(edge.arithId() == Arithmetics::SETNEG)
         {
            Datatypes::FieldTools::negative(*pOutVar);
         }

         if(edge.outId<int>() >= 0)
         {
            coord.communicator().storage<Dimensions::Transform::TRA2D>().holdBwd(*pOutVar, edge.outId<int>());
         } else
         {
            coord.communicator().transferBackward<Dimensions::Transform::TRA2D>(*pOutVar);
         }

      } else
      {
            coord.communicator().storage<Dimensions::Transform::TRA2D>().freeBwd(*pOutVar);
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD2D);
   }

}
}
