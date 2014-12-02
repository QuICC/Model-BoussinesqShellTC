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

namespace GeoMHDiSCC {

   namespace Transform {

      /**
       * @brief Struct to contain the full field transform steps
       *
       * The steps should respect the following naming structure:
       *    DO_???         : Step receives data, does computation, frees input and transfers output
       *    START_???      : Step receives data, does computation, holds input data and transfers output
       *    CONTINUE_???   : Step recovers data, does computations, holds input and transfers output
       *    FINISH_???     : Step recovers data, does computations, frees input and transfers output
       */
      namespace TransformSteps
      {
         template <Dimensions::Transform::Id TId> struct Forward;
         template <Dimensions::Transform::Id TId> struct Backward;

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
            static const Step STEP_VECTOR_ONE = NOTHING;

            /// Transform step for the second component of a vector
            static const Step STEP_VECTOR_TWO = NOTHING;

            /// Transform step for the third component of a vector
            static const Step STEP_VECTOR_THREE = NOTHING;


            /// First spectral field component used for vector transform
            static const FieldComponents::Spectral::Id SPECTOR_ONE = FieldComponents::Spectral::NOTUSED;

            /// Second spectral field component used for vector transform
            static const FieldComponents::Spectral::Id SPECTOR_TWO = FieldComponents::Spectral::NOTUSED;

            /// Third spectral field component used for vector transform
            static const FieldComponents::Spectral::Id SPECTOR_THREE = FieldComponents::Spectral::NOTUSED;

            /// Number of variable to transfer to next stage for scalar
            static const int SCALAR_VARIABLES = 1;

            /// Number of variable to transfer to next stage for vector
            static const int VECTOR_VARIABLES = 3;
         };

         /**
          * @brief Simple struct holding the backward transform steps for the second dimension
          */
         struct BackwardBase
         {
            /// Enum of the possible steps
            enum Step {NOTHING, DO_SCALAR, FINISH_SCALAR, CONTINUE_DIVR, FINISH_DIVR, START_GRAD, DO_GRAD};

            /// Transform step for a scalar
            static const Step STEP_SCALAR = DO_SCALAR;


            /// Transform step for the first component of a gradient
            static const Step STEP_GRAD_ONE = NOTHING;

            /// Transform step for the second component of a gradient
            static const Step STEP_GRAD_TWO = NOTHING;

            /// Transform step for the third component of a gradient
            static const Step STEP_GRAD_THREE = NOTHING;


            /// Transform step for the first component of a vector
            static const Step STEP_VECTOR_ONE = NOTHING;

            /// Transform step for the second component of a vector
            static const Step STEP_VECTOR_TWO = NOTHING;

            /// Transform step for the third component of a vector
            static const Step STEP_VECTOR_THREE = NOTHING;


            /// Transform step for the first component of a vector gradient
            static const Step STEP_VGRAD_ONE = NOTHING;

            /// Transform step for the second component of a vector gradient
            static const Step STEP_VGRAD_TWO = NOTHING;

            /// Transform step for the third component of a vector gradient
            static const Step STEP_VGRAD_THREE = NOTHING;


            /// Transform step for the first component of a curl
            static const Step STEP_CURL_ONE = NOTHING;

            /// Transform step for the second component of a curl
            static const Step STEP_CURL_TWO = NOTHING;

            /// Transform step for the third component of a curl
            static const Step STEP_CURL_THREE = NOTHING;


            /// First spectral field component used for vector
            static const FieldComponents::Spectral::Id SPECTOR_ONE = FieldComponents::Spectral::NOTUSED;

            /// Second spectral field component used for vector
            static const FieldComponents::Spectral::Id SPECTOR_TWO = FieldComponents::Spectral::NOTUSED;

            /// Third spectral field component used for vector
            static const FieldComponents::Spectral::Id SPECTOR_THREE = FieldComponents::Spectral::NOTUSED;


            /// First spectral field component used for curl
            static const FieldComponents::Spectral::Id SPECURL_ONE = FieldComponents::Spectral::NOTUSED;

            /// Second spectral field component used for curl
            static const FieldComponents::Spectral::Id SPECURL_TWO = FieldComponents::Spectral::NOTUSED;

            /// Third spectral field component used for curl
            static const FieldComponents::Spectral::Id SPECURL_THREE = FieldComponents::Spectral::NOTUSED;

            /// Number of variable to transfer to next stage for scalar
            static const int SCALAR_VARIABLES = 1;

            /// Number of variable to transfer to next stage for vector
            static const int VECTOR_VARIABLES = 3;

            /// Number of variable to transfer to next stage for curl
            static const int CURL_VARIABLES = 3;
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
         };

         /**
          * @brief Simple struct holding the backward transform steps for the first dimension
          */
         template<>  struct Backward<Dimensions::Transform::TRA1D>: public BackwardBase
         {
            /// Transform step for the first component of a gradient
            static const Step STEP_GRAD_ONE = START_GRAD;

            /// Transform step for the second component of a gradient
            static const Step STEP_GRAD_TWO = FINISH_SCALAR;

            /// Number of variable to transfer to next stage for gradient
            static const int GRAD_VARIABLES = 2;
         };

         /**
          * @brief Simple struct holding the backward transform steps for the second dimension
          */
         template<>  struct Backward<Dimensions::Transform::TRA2D>: public BackwardBase
         {
            /// Transform step for the first component of a gradient
            static const Step STEP_GRAD_ONE = DO_SCALAR;

            /// Transform step for the second component of a gradient
            static const Step STEP_GRAD_TWO = START_GRAD;

            /// Transform step for the third component of a gradient
            static const Step STEP_GRAD_THREE = FINISH_SCALAR;

            /// Number of variable to transfer to next stage for gradient
            static const int GRAD_VARIABLES = 3;
         };

         /**
          * @brief Simple struct holding the backward transform steps for the third dimension
          */
         template<>  struct Backward<Dimensions::Transform::TRA3D>: public BackwardBase
         {
            /// Transform step for the first component of a gradient
            static const Step STEP_GRAD_ONE = DO_SCALAR;

            /// Transform step for the second component of a gradient
            static const Step STEP_GRAD_TWO = DO_SCALAR;

            /// Transform step for the third component of a gradient
            static const Step STEP_GRAD_THREE = DO_GRAD;

            /// Number of variable to transfer to next stage for gradient
            static const int GRAD_VARIABLES = 3;
         };
      }
   }
}

#endif // CARTESIANTRANSFORMSTEPS_HPP
