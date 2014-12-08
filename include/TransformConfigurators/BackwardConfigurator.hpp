/** 
 * @file BackwardConfigurator.hpp
 * @brief This class defines the base operations for a backward transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BACKWARDCONFIGURATOR_HPP
#define BACKWARDCONFIGURATOR_HPP

// Configuration includes
// 
#include "Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TransformConfigurators/TransformStepsMacro.h"
#include "TransformConfigurators/ProjectorTree.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a backward transform
    */
   class BackwardConfigurator
   {
      public:

      protected:
         /**
          * @brief Prepare computation of projection for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void prepareSpectral(const ProjectorTree& tree, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of projection for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void prepareSpectral(const ProjectorTree& tree, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of physical projection for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void preparePhysical(const ProjectorTree& tree, const ProjectorTree::Projector3DEdge& edge, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of physical projection for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void preparePhysical(const ProjectorTree& tree, const ProjectorTree::Projector3DEdge& edge, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Compute the projection transform of the first dimension
          *
          * @param edge    Transform tree edge
          * @param coord   Transform coordinator
          * @param hold    Hold input data?
          */
         static void project1D(const ProjectorTree::Projector1DEdge& edge, TransformCoordinatorType& coord, const bool hold);

         /**
          * @brief Compute the projection transform of the second dimension
          *
          * @param edge    Transform tree edge
          * @param coord   Transform coordinator
          * @param recover Recover input data?
          * @param hold    Hold input data?
          */
         static void project2D(const ProjectorTree::Projector2DEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold);

         /**
          * @brief Compute the projection transform of the third dimension
          *
          * @param edge    Transform tree edge
          * @param coord   Transform coordinator
          * @param recover Recover input data?
          * @param hold    Hold input data?
          */
         static void project3D(const ProjectorTree::Projector3DEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold);

         /**
          * @brief Empty constructor
          */
         BackwardConfigurator();

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardConfigurator();

      private: 
   };

}
}

#endif // BACKWARDCONFIGURATOR_HPP
