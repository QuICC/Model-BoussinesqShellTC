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
#include "ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Transform {

   void BackwardConfigurator2D::prepareSpectral(const TransformTree& tree, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Safety assert
      assert(tree.comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);
      DetailedProfilerMacro_start(ProfilerMacro::BWDDEALIAS);

      // Put scalar into temporary hold storage
      coord.communicator().dealiasSpectral(rScalar.rDom(0).rTotal());

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDDEALIAS);
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   void BackwardConfigurator2D::prepareSpectral(const TransformTree& tree, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Safety assert
      assert(tree.comp<FieldComponents::Spectral::Id>() != FieldComponents::Spectral::SCALAR);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Put scalar into temporary hold storage
      coord.communicator().dealiasSpectral(rVector.rDom(0).rTotal().rComp(tree.comp<FieldComponents::Spectral::Id>()));

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   void BackwardConfigurator2D::preparePhysical(const TransformTree& tree, const TransformTreeEdge& edge, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
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
         coord.communicator().holdPhysical(rScalar.rDom(0).rGrad().rComp(edge.outId<FieldComponents::Physical::Id>()));

      // Put 2nd order gradient component into temporary hold storage
      } else if(edge.fieldId() == FieldType::GRADIENT2)
      {
         coord.communicator().holdPhysical(rScalar.rDom(0).rGrad2().rComp(edge.outId<FieldComponents::Physical::Id>(0),edge.outId<FieldComponents::Physical::Id>(1)));
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDND);
   }

   void BackwardConfigurator2D::preparePhysical(const TransformTree& tree, const TransformTreeEdge& edge, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDND);

      // Put vector component into temporary hold storage
      if(edge.fieldId() == FieldType::VECTOR)
      {
         coord.communicator().holdPhysical(rVector.rDom(0).rPhys().rComp(edge.outId<FieldComponents::Physical::Id>()));

      // Put vector gradient component into temporary hold storage
      } else if(edge.fieldId() == FieldType::GRADIENT)
      {
         coord.communicator().holdPhysical(rVector.rDom(0).rGrad(tree.comp<FieldComponents::Spectral::Id>()).rComp(edge.outId<FieldComponents::Physical::Id>()));

      // Put curl component into temporary hold storage
      } else if(edge.fieldId() == FieldType::CURL)
      {
         coord.communicator().holdPhysical(rVector.rDom(0).rCurl().rComp(edge.outId<FieldComponents::Physical::Id>()));
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDND);
   }

   void BackwardConfigurator2D::project1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Project 1D", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get output storage
      TransformCoordinatorType::CommunicatorType::Fwd1DType *pOutVar = 0;
      TransformCoordinatorType::CommunicatorType::Fwd1DType *pRecOutVar = 0;
      if(edge.recoverOutId() >= 0)
      {
         pOutVar = &coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();
         pRecOutVar = &coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverFwd(edge.recoverOutId());
      } else
      {
         pOutVar = &coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();
      }

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1DTRA);

      // Compute projection transform for first dimension 
      coord.transform1D().project(pOutVar->rData(), rInVar.data(), edge.opId<TransformSelector<Dimensions::Transform::TRA1D>::Type::ProjectorType::Id>());
      
      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1DTRA);

      // Hold spectral input
      if(edge.holdInput())
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().holdBwd(rInVar);

      // Free spectral input
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);
      }

      // Combine recovered output with new calculation
      if(pRecOutVar != 0)
      {
         Datatypes::FieldTools::combine(*pRecOutVar, *pOutVar, edge.combinedArithId());

         if(edge.combinedOutId() >= 0)
         {
            coord.communicator().storage<Dimensions::Transform::TRA1D>().holdFwd(*pRecOutVar, edge.combinedOutId());
         } else
         {
            coord.communicator().transferForward<Dimensions::Transform::TRA1D>(*pRecOutVar);
         }

      // Hold data for combination   
      } else if(edge.combinedOutId() >= 0)
      {
         if(edge.combinedArithId() == Arithmetics::SETNEG)
         {
            Datatypes::FieldTools::negative(*pOutVar);
         }

         coord.communicator().storage<Dimensions::Transform::TRA1D>().holdFwd(*pOutVar, edge.combinedOutId());
      }

      // Transfer calculation
      if(edge.arithId() != Arithmetics::NONE)
      {
         if(edge.arithId() == Arithmetics::SETNEG)
         {
            Datatypes::FieldTools::negative(*pOutVar);
         }

         coord.communicator().transferForward<Dimensions::Transform::TRA1D>(*pOutVar);

      } else if(edge.combinedOutId() < 0 || pRecOutVar != 0)
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(*pOutVar);
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   void BackwardConfigurator2D::projectND(const TransformTreeEdge& edge, TransformCoordinatorType& coord)
   {
      // Debugger message
      DebuggerMacro_msg("Project ND", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDND);

      TransformCoordinatorType::CommunicatorType::BwdNDType *pInVar;

      // Get the input data from hold
      if(edge.recoverInput())
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRAND>().recoverBwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveBackward<Dimensions::Transform::TRAND>();
      }

      // Get output storage
      TransformCoordinatorType::CommunicatorType::FwdNDType *pOutVar = 0;
      TransformCoordinatorType::CommunicatorType::FwdNDType *pRecOutVar = 0;
      if(edge.arithId() == Arithmetics::SET || edge.arithId() == Arithmetics::SETNEG)
      {
         pOutVar = &coord.communicator().storage<Dimensions::Transform::TRAND>().recoverFwd();
      } else
      {
         pOutVar = &coord.communicator().storage<Dimensions::Transform::TRAND>().provideFwd();
         pRecOutVar = &coord.communicator().storage<Dimensions::Transform::TRAND>().recoverFwd();
      }

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDNDTRA);

      // Compute projection transform for third dimension 
      coord.transformND().project(pOutVar->rData(), pInVar->data(), edge.opId<TransformSelector<Dimensions::Transform::TRAND>::Type::ProjectorType::Id>());

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDNDTRA);

      // Hold temporary storage
      if(edge.holdInput())
      {
         coord.communicator().storage<Dimensions::Transform::TRAND>().holdBwd(*pInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRAND>().freeBwd(*pInVar);
      }

      // Combine recovered output with new calculation
      if(pRecOutVar != 0)
      {
         Datatypes::FieldTools::combine(*pRecOutVar, *pOutVar, edge.arithId());

         coord.communicator().transferForward<Dimensions::Transform::TRAND>(*pRecOutVar);
         coord.communicator().storage<Dimensions::Transform::TRAND>().freeFwd(*pOutVar);
      } else
      {
         if(edge.arithId() == Arithmetics::SETNEG)
         {
            Datatypes::FieldTools::negative(*pOutVar);
         }

         coord.communicator().transferForward<Dimensions::Transform::TRAND>(*pOutVar);
      }

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDND);
   }

}
}
