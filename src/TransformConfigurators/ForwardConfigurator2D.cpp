/** 
 * @file ForwardConfigurator2D.cpp
 * @brief Source of the implementation of the base forward configurator in 2D space
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
#include "TransformConfigurators/ForwardConfigurator2D.hpp"

// Project includes
//
#include "ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Transform {

   void ForwardConfigurator2D::nonlinearTerm(const TransformTree& tree, Equations::SharedIEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::NONLINEAR);

      // Get physical storage
      TransformCoordinatorType::CommunicatorType::FwdNDType &rNLComp = coord.communicator().providePhysical();

      // Compute nonlinear term component
      spEquation->computeNonlinear(rNLComp, tree.comp<FieldComponents::Physical::Id>());
      spEquation->useNonlinear(rNLComp, tree.comp<FieldComponents::Physical::Id>());

      // Transfer physical storage to next step
      coord.communicator().holdPhysical(rNLComp);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::NONLINEAR);
   }

   void ForwardConfigurator2D::integrateND(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate ND", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWDND);

      // Get recover input data from hold
      TransformCoordinatorType::CommunicatorType::FwdNDType &rInVar = coord.communicator().receiveForward<Dimensions::Transform::TRAND>();

      // Get output storage
      TransformCoordinatorType::CommunicatorType::BwdNDType *pOutVar = 0;
      TransformCoordinatorType::CommunicatorType::BwdNDType *pRecOutVar = 0;
      if(edge.recoverOutId() >= 0)
      {
         pOutVar = &coord.communicator().storage<Dimensions::Transform::TRAND>().provideBwd();
         pRecOutVar = &coord.communicator().storage<Dimensions::Transform::TRAND>().recoverBwd(edge.recoverOutId());
      } else
      {
         pOutVar = &coord.communicator().storage<Dimensions::Transform::TRAND>().provideBwd();
      }

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWDNDTRA);

      // Compute integration transform of third dimension
      coord.transformND().integrate(pOutVar->rData(), rInVar.data(), edge.opId<TransformSelector<Dimensions::Transform::TRAND>::Type::IntegratorType::Id>());

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDNDTRA);

      // Hold temporary input storage
      if(edge.holdInput())
      {
         coord.communicator().storage<Dimensions::Transform::TRAND>().holdFwd(rInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRAND>().freeFwd(rInVar);
      }

      // Combine recovered output with new calculation
      if(pRecOutVar != 0)
      {
         Datatypes::FieldTools::combine(*pRecOutVar, *pOutVar, edge.combinedArithId());

         if(edge.combinedOutId() >= 0)
         {
            coord.communicator().storage<Dimensions::Transform::TRAND>().holdBwd(*pRecOutVar, edge.combinedOutId());
         } else
         {
            coord.communicator().transferBackward<Dimensions::Transform::TRAND>(*pRecOutVar);
         }

      // Hold data for combination   
      } else if(edge.combinedOutId() >= 0)
      {
         if(edge.combinedArithId() == Arithmetics::SETNEG)
         {
            Datatypes::FieldTools::negative(*pOutVar);
         }

         coord.communicator().storage<Dimensions::Transform::TRAND>().holdBwd(*pOutVar, edge.combinedOutId());
      }

      // Transfer calculation
      if(edge.arithId() != Arithmetics::NONE)
      {
         if(edge.arithId() == Arithmetics::SETNEG)
         {
            Datatypes::FieldTools::negative(*pOutVar);
         }

         coord.communicator().transferBackward<Dimensions::Transform::TRAND>(*pOutVar);

      } else if(edge.combinedOutId() < 0 || pRecOutVar != 0)
      {
         coord.communicator().storage<Dimensions::Transform::TRAND>().freeBwd(*pOutVar);
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDND);
   }

   void ForwardConfigurator2D::integrate1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 1D", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD1D);

      TransformCoordinatorType::CommunicatorType::Fwd1DType* pInVar;

      // Recover hold input data
      if(edge.recoverInput())
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverFwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveForward<Dimensions::Transform::TRA1D>();
      }

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideBwd();

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD1DTRA);

      // Compute integration transform of first dimension
      coord.transform1D().integrate(rOutVar.rData(), pInVar->data(), edge.opId<TransformSelector<Dimensions::Transform::TRA1D>::Type::IntegratorType::Id>());

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD1DTRA);

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().holdFwd(*pInVar);

      // Free temporary storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(*pInVar);
      }

      // Hold temporary storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().holdBwd(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD1D);
   }

}
}
