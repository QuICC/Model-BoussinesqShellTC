/** 
 * @file BackwardConfigurator2D.cpp
 * @brief Source of the implementation of the base backward configurator
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
#include "TransformConfigurators/BackwardConfigurator2D.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void BackwardConfigurator2D::prepareSpectral(const ProjectorTree& tree, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Safety assert
      assert(tree.comp() == FieldComponents::Spectral::SCALAR);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Put scalar into temporary hold storage
      coord.communicator().dealiasSpectral(rScalar.rDom(0).rTotal());

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   void BackwardConfigurator2D::prepareSpectral(const ProjectorTree& tree, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Safety assert
      assert(tree.comp() != FieldComponents::Spectral::SCALAR);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Put scalar into temporary hold storage
      coord.communicator().dealiasSpectral(rVector.rDom(0).rTotal().rComp(tree.comp()));

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   void BackwardConfigurator2D::preparePhysical(const ProjectorTree& tree, const ProjectorPhysEdge& edge, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDND);

      // Put scalar into temporary hold storage
      if(edge.fieldId() == FieldType::SCALAR)
      {
         coord.communicator().holdPhysical(rScalar.rDom(0).rPhys());

      // Put gradient component into temporary hold storage
      } else if(edge.fieldId() == FieldType::GRADIENT)
      {
         coord.communicator().holdPhysical(rScalar.rDom(0).rGrad().rComp(edge.physId()));
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDND);
   }

   void BackwardConfigurator2D::preparePhysical(const ProjectorTree& tree, const ProjectorPhysEdge& edge, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDND);

      // Put vector component into temporary hold storage
      if(edge.fieldId() == FieldType::VECTOR)
      {
         coord.communicator().holdPhysical(rVector.rDom(0).rPhys().rComp(edge.physId()));

      // Put vector gradient component into temporary hold storage
      } else if(edge.fieldId() == FieldType::GRADIENT)
      {
         coord.communicator().holdPhysical(rVector.rDom(0).rGrad(tree.comp()).rComp(edge.physId()));

      // Put curl component into temporary hold storage
      } else if(edge.fieldId() == FieldType::CURL)
      {
         coord.communicator().holdPhysical(rVector.rDom(0).rCurl().rComp(edge.physId()));
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDND);
   }

   void BackwardConfigurator2D::project1D(const ProjectorSpecEdge& edge, TransformCoordinatorType& coord, const bool hold)
   {
      // Debugger message
      DebuggerMacro_msg("Project 1D", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRA);

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVar.rData(), rInVar.data(), edge.opId(), Arithmetics::SET);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRA);

      // Hold spectral input
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().holdBwd(rInVar);

      // Free spectral input
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);
      }

      // Transfer output data to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA1D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   void BackwardConfigurator2D::projectND(const ProjectorPhysEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold)
   {
      // Debugger message
      DebuggerMacro_msg("Project ND", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDND);

      TransformCoordinatorType::CommunicatorType::BwdNDType *pInVar;

      // Get the input data from hold
      if(recover)
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRAND>().recoverBwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveBackward<Dimensions::Transform::TRAND>();
      }

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::FwdNDType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRAND>().recoverFwd();

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDNDTRA);

      // Compute projection transform for third dimension 
      coord.transformND().project(rOutVar.rData(), pInVar->data(), edge.opId(),  edge.arithId());

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDNDTRA);

      // Hold temporary storage
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRAND>().holdBwd(*pInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRAND>().freeBwd(*pInVar);
      }

      // Hold temporary storage
      coord.communicator().transferForward<Dimensions::Transform::TRAND>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDND);
   }

}
}
