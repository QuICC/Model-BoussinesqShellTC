/** 
 * @file BackwardConfigurator2D.hpp
 * @brief This class defines the base operations for a backward transform in 2D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BACKWARDCONFIGURATOR2D_HPP
#define BACKWARDCONFIGURATOR2D_HPP

// Configuration includes
// 
#include "Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/TransformCommSelector.hpp"
#include "TransformConfigurators/TransformTree.hpp"
#include "TransformConfigurators/TransformStepsMacro.h"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a backward transform in 2D space
    */
   class BackwardConfigurator2D
   {
      public:

      protected:
         /**
          * @brief Prepare computation of projection for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void prepareSpectral(const TransformTree& tree, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of projection for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void prepareSpectral(const TransformTree& tree, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of physical projection for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void preparePhysical(const TransformTree& tree, const TransformTreeEdge& edge, Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of physical projection for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void preparePhysical(const TransformTree& tree, const TransformTreeEdge& edge, Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Compute the projection transform of the first dimension
          *
          * @param edge    Transform tree edge
          * @param coord   Transform coordinator
          */
         static void project1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Compute the projection transform of the third dimension
          *
          * @param edge    Transform tree edge
          * @param coord   Transform coordinator
          */
         static void projectND(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         BackwardConfigurator2D();

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardConfigurator2D();

      private: 
   };

#ifdef GEOMHDISCC_SPATIALDIMENSION_2D
   typedef BackwardConfigurator2D BackwardConfigurator;
#endif //GEOMHDISCC_SPATIALDIMENSION_2D

}
}

#endif // BACKWARDCONFIGURATOR2D_HPP
