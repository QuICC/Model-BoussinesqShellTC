/** 
 * @file BackwardConfigurator.cpp
 * @brief Source of the implementation of the base backward configurator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/BackwardConfigurator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void BackwardConfigurator::prepareSpectral(const ProjectorTree& tree, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
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

   void BackwardConfigurator::prepareSpectral(const ProjectorTree& tree, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
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

   void BackwardConfigurator::preparePhysical(const ProjectorTree& tree, const ProjectorTree::Projector3DEdge& edge, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

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
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

   void BackwardConfigurator::preparePhysical(const ProjectorTree& tree, const ProjectorTree::Projector3DEdge& edge, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

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
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

   void BackwardConfigurator::project1D(const ProjectorTree::Projector1DEdge& edge, TransformCoordinatorType& coord, const bool hold)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), edge.opId());

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

   void BackwardConfigurator::project2D(const ProjectorTree::Projector2DEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd2DType* pInVar;
      if(recover)
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().recoverBwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveBackward<Dimensions::Transform::TRA2D>();
      }

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd2DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA2D>().provideFwd();

      // Compute projection transform for second dimension 
      coord.transform2D().project<Arithmetics::SET>(rOutVar.rData(), pInVar->data(), edge.opId());

      // Hold temporary storage
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().holdBwd(*pInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().freeBwd(*pInVar);
      }

      // Transfer output data to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA2D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2D);
   }

   void BackwardConfigurator::project3D(const ProjectorTree::Projector3DEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

      TransformCoordinatorType::CommunicatorType::Bwd3DType *pInVar;

      // Get the input data from hold
      if(recover)
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRA3D>().recoverBwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveBackward<Dimensions::Transform::TRA3D>();
      }

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA3D>().recoverFwd();

      // Compute projection transform for third dimension 
      coord.transform3D().project<Arithmetics::SET>(rOutVar.rData(), pInVar->data(), edge.opId());

      // Hold temporary storage
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRA3D>().holdBwd(*pInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA3D>().freeBwd(*pInVar);
      }

      // Hold temporary storage
      coord.communicator().transferForward<Dimensions::Transform::TRA3D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

}
}
