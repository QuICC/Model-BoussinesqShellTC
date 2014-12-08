/** 
 * @file ForwardSerialConfigurator.cpp
 * @brief Source of the implementation of the forward serial configurator
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
#include "TransformConfigurators/ForwardSerialConfigurator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void ForwardSerialConfigurator::firstStep(const IntegratorTree& tree, Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      IntegratorTree::Integrator1DEdge_iterator it1D;
      IntegratorTree::Integrator2DEdge_iterator it2D;
      IntegratorTree::Integrator3DEdge_iterator it3D;

      // Ranges for the vector of edges for the three transforms
      IntegratorTree::Integrator1DEdge_range range1D;
      IntegratorTree::Integrator2DEdge_range range2D;
      IntegratorTree::Integrator3DEdge_range range3D = tree.edgeRange();

      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm(tree, spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Loop over first transform
      int hold3D = std::distance(range3D.first, range3D.second) - 1;
      for(it3D = range3D.first; it3D != range3D.second; ++it3D, --hold3D)
      {
         // Compute third transform
         ForwardConfigurator::integrate3D(*it3D, coord, hold3D);

         range2D = it3D->edgeRange();
         int recover2D = 0;
         int hold2D = std::distance(range2D.first, range2D.second) - 1;
         for(it2D = range2D.first; it2D != range2D.second; ++it2D, ++recover2D, --hold2D)
         {
            // Compute second transform
            ForwardConfigurator::integrate2D(*it2D, coord, recover2D, hold2D);

            range1D = it2D->edgeRange();
            int recover1D = 0;
            int hold1D = std::distance(range1D.first, range1D.second) - 1;
            for(it1D = range1D.first; it1D != range1D.second; ++it1D, ++recover1D, --hold1D)
            {
               // Compute third transform
               ForwardConfigurator::integrate1D(*it1D, coord, recover1D, hold1D);
            }
         }

         // Stop profiler
         ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
      }
   }

   void ForwardSerialConfigurator::firstStep(const IntegratorTree& tree, Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
   }

   void ForwardSerialConfigurator::secondStep(const IntegratorTree& tree, Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a second step
   }

   void ForwardSerialConfigurator::secondStep(const IntegratorTree& tree, Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a second step
   }
   
   void ForwardSerialConfigurator::lastStep(const IntegratorTree& tree, Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a last step
   }

   void ForwardSerialConfigurator::lastStep(const IntegratorTree& tree, Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a last step
   }

}
}
