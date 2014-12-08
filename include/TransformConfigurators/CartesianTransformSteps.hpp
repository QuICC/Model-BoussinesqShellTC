/** 
 * @file CartesianTransformSteps.hpp
 * @brief Definition of some useful enums for steps involved in a full cartesian transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIANTRANSFORMSTEPS_HPP
#define CARTESIANTRANSFORMSTEPS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "TransformConfigurators/ProjectorBranch.hpp"
#include "TransformConfigurators/IntegratorBranch.hpp"

namespace GeoMHDiSCC {

   namespace Transform {

      /**
       * @brief Struct to contain the full field transform steps
       */
      namespace TransformSteps
      {
         /**
          * @brief Generate the list of branches in scalar integration transform
          */
         std::vector<IntegratorBranch>  forwardScalar();

         /**
          * @brief Generate the list of branches in vector integration transform
          */
         std::vector<IntegratorBranch>  forwardVector();

         /**
          * @brief Generate the list of branches in scalar projection transform
          */
         std::vector<ProjectorBranch>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in scalar gradient transform
          */
         std::vector<ProjectorBranch>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in vector projection transform
          */
         std::vector<ProjectorBranch>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in vector gradient transform
          */
         std::vector<ProjectorBranch>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in vector curl transform
          */
         std::vector<ProjectorBranch>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in vector divergence transform
          */
         std::vector<ProjectorBranch>  backwardDivergence();

         template <Dimensions::Transform::Id TId> struct Forward;

         /**
          * @brief Simple struct holding the forward transform steps for the first dimension
          */
         struct ForwardBase
         {
            /// Enum of the possible steps
            enum Step {NOTHING, DO_SCALAR};

            /// Transform step for a scalar
            static const Step STEP_SCALAR = DO_SCALAR;


            /// Transform step for the first component of a vector
            static const Step STEP_VECTOR_ONE = DO_SCALAR;

            /// Transform step for the second component of a vector
            static const Step STEP_VECTOR_TWO = DO_SCALAR;

            /// Transform step for the third component of a vector
            static const Step STEP_VECTOR_THREE = DO_SCALAR;


            /// First spectral field component used for vector transform
            static const FieldComponents::Spectral::Id SPECTOR_ONE = FieldComponents::Spectral::ONE;

            /// Second spectral field component used for vector transform
            static const FieldComponents::Spectral::Id SPECTOR_TWO = FieldComponents::Spectral::TWO;

            /// Third spectral field component used for vector transform
            static const FieldComponents::Spectral::Id SPECTOR_THREE = FieldComponents::Spectral::THREE;


            /// Number of variable to transfer to next stage for scalar
            static const int SCALAR_VARIABLES = 1;

            /// Number of variable to transfer to next stage for vector
            static const int VECTOR_VARIABLES = 3;
         };

         /**
          * @brief Simple struct holding the forward transform steps for the first dimension
          */
         template<>  struct Forward<Dimensions::Transform::TRA1D>: public ForwardBase
         {
         };

         /**
          * @brief Simple struct holding the forward transform steps for the second dimension
          */
         template<>  struct Forward<Dimensions::Transform::TRA2D>: public ForwardBase
         {
         };

         /**
          * @brief Simple struct holding the forward transform steps for the third dimension
          */
         template<>  struct Forward<Dimensions::Transform::TRA3D>: public ForwardBase
         {
         };

         /**
          * @brief Simple struct holding the physical components order
          */
         struct PhysicalBase
         {
            /// First physical field component used for gradient
            static const FieldComponents::Physical::Id GRAD_ONE = FieldComponents::Physical::NOTUSED;

            /// Second physical field component used for gradient
            static const FieldComponents::Physical::Id GRAD_TWO = FieldComponents::Physical::NOTUSED;

            /// Third physical field component used for gradient
            static const FieldComponents::Physical::Id GRAD_THREE = FieldComponents::Physical::NOTUSED;


            /// First physical field component used for vector
            static const FieldComponents::Physical::Id NONLINEAR_ONE = FieldComponents::Physical::NOTUSED;

            /// Second physical field component used for vector
            static const FieldComponents::Physical::Id NONLINEAR_TWO = FieldComponents::Physical::NOTUSED;

            /// Third physical field component used for vector
            static const FieldComponents::Physical::Id NONLINEAR_THREE = FieldComponents::Physical::NOTUSED;


            /// First physical field component used for vector
            static const FieldComponents::Physical::Id VECTOR_ONE = FieldComponents::Physical::NOTUSED;

            /// Second physical field component used for vector
            static const FieldComponents::Physical::Id VECTOR_TWO = FieldComponents::Physical::NOTUSED;

            /// Third physical field component used for vector
            static const FieldComponents::Physical::Id VECTOR_THREE = FieldComponents::Physical::NOTUSED;


            /// First physical field component used for vector gradient
            static const FieldComponents::Physical::Id VGRAD_ONE = FieldComponents::Physical::NOTUSED;

            /// Second physical field component used for vector gradient
            static const FieldComponents::Physical::Id VGRAD_TWO = FieldComponents::Physical::NOTUSED;

            /// Third physical field component used for vector gradient
            static const FieldComponents::Physical::Id VGRAD_THREE = FieldComponents::Physical::NOTUSED;


            /// First physical field component used for curl
            static const FieldComponents::Physical::Id CURL_ONE = FieldComponents::Physical::NOTUSED;

            /// Second physical field component used for curl
            static const FieldComponents::Physical::Id CURL_TWO = FieldComponents::Physical::NOTUSED;

            /// Third physical field component used for curl
            static const FieldComponents::Physical::Id CURL_THREE = FieldComponents::Physical::NOTUSED;
         };

         /**
          * @brief Simple struct holding the physical components order
          */
         struct Physical: public PhysicalBase
         {
            /// First physical field component used for gradient
            static const FieldComponents::Physical::Id GRAD_ONE = FieldComponents::Physical::ONE;

            /// Second physical field component used for gradient
            static const FieldComponents::Physical::Id GRAD_TWO = FieldComponents::Physical::TWO;

            /// Third physical field component used for gradient
            static const FieldComponents::Physical::Id GRAD_THREE = FieldComponents::Physical::THREE;


            /// First physical field component used for vector
            static const FieldComponents::Physical::Id VECTOR_ONE = FieldComponents::Physical::ONE;

            /// Second physical field component used for vector
            static const FieldComponents::Physical::Id VECTOR_TWO = FieldComponents::Physical::TWO;

            /// Third physical field component used for vector
            static const FieldComponents::Physical::Id VECTOR_THREE = FieldComponents::Physical::THREE;
         };
      }
   }
}

#endif // CARTESIANTRANSFORMSTEPS_HPP
