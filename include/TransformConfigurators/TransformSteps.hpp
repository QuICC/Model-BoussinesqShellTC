/** \file TransformSteps.hpp
 *  \brief Definition of some useful enums for steps involved in a full transform
 */

#ifndef TRANSFORMSTEPS_HPP
#define TRANSFORMSTEPS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/FieldComponents.hpp"

namespace GeoMHDiSCC {

   namespace Transform {

      /**
       * @brief Struct to contain the full field transform steps
       */
      namespace TransformSteps
      {
         template <Dimensions::Transform::Id TId> struct Forward;
         template <Dimensions::Transform::Id TId> struct Backward;

         /**
          * @brief Simple struct holding the forward transform steps for the first dimension
          */
         template<>  struct Forward<Dimensions::Transform::TRA1D>
         {
            /// Enum of the possible steps
            enum Step {NOTHING, PROJECTSCALAR};

            /// Transform step for a scalar
            static const Step SCALAR = NOTHING;


            /// Transform step for the first component of a vector
            static const Step VECTORONE = NOTHING;

            /// Transform step for the second component of a vector
            static const Step VECTORTWO = NOTHING;

            /// Transform step for the third component of a vector
            static const Step VECTORTHREE = NOTHING;

            /// First spectral field component used for vector
            static const FieldComponents::Spectral::Id SPECVECT1D = FieldComponents::Spectral::NONE;

            /// Second spectral field component used for vector
            static const FieldComponents::Spectral::Id SPECVECT2D = FieldComponents::Spectral::NONE;

            /// Third spectral field component used for vector
            static const FieldComponents::Spectral::Id SPECVECT3D = FieldComponents::Spectral::NONE;
         };

         /**
          * @brief Simple struct holding the forward transform steps for the second dimension
          */
         template<>  struct Forward<Dimensions::Transform::TRA2D>
         {
            /// Enum of the possible steps
            enum Step {NOTHING, PROJECTSCALAR};

            /// Transform step for a scalar
            static const Step SCALAR = SCALAR;


            /// Transform step for the first component of a vector
            static const TransformSteps::Id::Fwd2D::Step VECTOR1D = TransformSteps::Id::Fwd2D::NOTHING;

            /// Transform step for the second component of a vector
            static const TransformSteps::Id::Fwd2D::Step VECTOR2D = TransformSteps::Id::Fwd2D::NOTHING;

            /// Transform step for the third component of a vector
            static const TransformSteps::Id::Fwd2D::Step VECTOR3D = TransformSteps::Id::Fwd2D::NOTHING;
         };

         /**
          * @brief Simple struct holding the forward transform steps for the third dimension
          */
         template<>  struct Forward<Dimensions::Transform::TRA3D>
         {
            /// Enum of the possible steps
            enum Step {NOTHING, SCALAR};

            /// Transform step for a scalar
            static const Step SCALAR = SCALAR;


            /// Transform step for the first component of a vector
            static const TransformSteps::Id::Fwd3D::Step VECTOR1D = TransformSteps::Id::Fwd3D::NOTHING;

            /// Transform step for the second component of a vector
            static const TransformSteps::Id::Fwd3D::Step VECTOR2D = TransformSteps::Id::Fwd3D::NOTHING;

            /// Transform step for the third component of a vector
            static const TransformSteps::Id::Fwd3D::Step VECTOR3D = TransformSteps::Id::Fwd3D::NOTHING;
         };

         /**
          * @brief Simple struct holding the physical components order
          */
         struct Physical {
            /// First physical field component used for gradient
            static const FieldComponents::Physical::Id GRAD1D = FieldComponents::Physical::ONE;

            /// Second physical field component used for gradient
            static const FieldComponents::Physical::Id GRAD2D = FieldComponents::Physical::TWO;

            /// Third physical field component used for gradient
            static const FieldComponents::Physical::Id GRAD3D = FieldComponents::Physical::THREE;


            /// First physical field component used for vector
            static const FieldComponents::Physical::Id NONLINEAR1D = FieldComponents::Physical::TWO;

            /// Second physical field component used for vector
            static const FieldComponents::Physical::Id NONLINEAR2D = FieldComponents::Physical::THREE;

            /// Third physical field component used for vector
            static const FieldComponents::Physical::Id NONLINEAR3D = FieldComponents::Physical::ONE;


            /// First physical field component used for vector
            static const FieldComponents::Physical::Id VECTOR1D = FieldComponents::Physical::ONE;

            /// Second physical field component used for vector
            static const FieldComponents::Physical::Id VECTOR2D = FieldComponents::Physical::TWO;

            /// Third physical field component used for vector
            static const FieldComponents::Physical::Id VECTOR3D = FieldComponents::Physical::ONE;


            /// First physical field component used for curl
            static const FieldComponents::Physical::Id CURL1D = FieldComponents::Physical::ONE;

            /// Second physical field component used for curl
            static const FieldComponents::Physical::Id CURL2D = FieldComponents::Physical::TWO;

            /// Third physical field component used for curl
            static const FieldComponents::Physical::Id CURL3D = FieldComponents::Physical::THREE;
         };

         /**
          * @brief Simple struct holding the backward transform steps for the first dimension
          */
         template<>  struct Backward<Dimensions::Transform::TRA1D>
         {
            /// Enum of the possible steps
            enum Step {NONE, SCALAR};

            /// Transform step for a scalar
            static const Step SCALAR = SCALAR;


            /// Transform step for the first component of a gradient
            static const Step GRAD1D = TransformSteps::Id::Bwd1D::GRADX;

            /// Transform step for the second component of a gradient
            static const Step GRAD2D = TransformSteps::Id::Bwd1D::SCALAR;

            /// Transform step for the third component of a gradient
            static const Step GRAD3D = TransformSteps::Id::Bwd1D::NONE;


            /// Transform step for the first component of a vector
            static const Step VECTOR1D = TransformSteps::Id::Bwd1D::NONE;

            /// Transform step for the second component of a vector
            static const Step VECTOR2D = TransformSteps::Id::Bwd1D::NONE;

            /// Transform step for the third component of a vector
            static const Step VECTOR3D = TransformSteps::Id::Bwd1D::NONE;

            /// First spectral field component used for vector
            static const FieldComponents::Spectral::Id SPECVECT1D = FieldComponents::Spectral::NONE;

            /// Second spectral field component used for vector
            static const FieldComponents::Spectral::Id SPECVECT2D = FieldComponents::Spectral::NONE;

            /// Third spectral field component used for vector
            static const FieldComponents::Spectral::Id SPECVECT3D = FieldComponents::Spectral::NONE;


            /// Transform step for the first component of a curl
            static const TransformSteps::Id::Bwd1D::Step CURL1D = TransformSteps::Id::Bwd1D::NONE;

            /// Transform step for the second component of a curl
            static const TransformSteps::Id::Bwd1D::Step CURL2D = TransformSteps::Id::Bwd1D::NONE;

            /// Transform step for the third component of a curl
            static const TransformSteps::Id::Bwd1D::Step CURL3D = TransformSteps::Id::Bwd1D::NONE;

            /// First spectral field component used for curl
            static const FieldComponents::Spectral::Id SPECCURL1D = FieldComponents::Spectral::NONE;

            /// Second spectral field component used for curl
            static const FieldComponents::Spectral::Id SPECCURL2D = FieldComponents::Spectral::NONE;

            /// Third spectral field component used for curl
            static const FieldComponents::Spectral::Id SPECCURL3D = FieldComponents::Spectral::NONE;
         };

         /**
          * @brief Simple struct holding the backward transform steps for the second dimension
          */
         template<>  struct Backward<Dimensions::Transform::TRA2D>
         {
            /// Enum of the possible steps
            enum Step {NONE, SCALAR};

            /// Transform step for a scalar
            static const Step SCALAR = TransformSteps::Id::Bwd2D::SCALAR;


            /// Transform step for the first component of a gradient
            static const Step GRAD1D = TransformSteps::Id::Bwd2D::SCALAR;

            /// Transform step for the second component of a gradient
            static const Step GRAD2D = TransformSteps::Id::Bwd2D::GRADY;

            /// Transform step for the third component of a gradient
            static const Step GRAD3D = TransformSteps::Id::Bwd2D::SCALAR;


            /// Transform step for the first component of a vector
            static const Step VECTOR1D = TransformSteps::Id::Bwd2D::NONE;

            /// Transform step for the second component of a vector
            static const Step VECTOR2D = TransformSteps::Id::Bwd2D::NONE;

            /// Transform step for the third component of a vector
            static const Step VECTOR3D = TransformSteps::Id::Bwd2D::NONE;


            /// Transform step for the first component of a curl
            static const Step CURL1D = TransformSteps::Id::Bwd2D::NONE;

            /// Transform step for the second component of a curl
            static const Step CURL2D = TransformSteps::Id::Bwd2D::NONE;

            /// Transform step for the third component of a curl
            static const Step CURL3D = TransformSteps::Id::Bwd2D::NONE;
         };

         /**
          * @brief Simple struct holding the backward transform steps for the third dimension
          */
         template<>  struct Backward<Dimensions::Transform::TRA3D>
         {
            /// Enum of the possible steps
            enum Step {NONE, SCALAR};

            /// Transform step for a scalar
            static const Step SCALAR = TransformSteps::Id::Bwd3D::SCALAR;


            /// Transform step for the first component of a gradient
            static const Step GRAD1D = TransformSteps::Id::Bwd3D::SCALAR;

            /// Transform step for the second component of a gradient
            static const Step GRAD2D = TransformSteps::Id::Bwd3D::SCALAR;

            /// Transform step for the third component of a gradient
            static const Step GRAD3D = TransformSteps::Id::Bwd3D::GRADZ;


            /// Transform step for the first component of a vector
            static const Step VECTOR1D = TransformSteps::Id::Bwd3D::NONE;

            /// Transform step for the second component of a vector
            static const Step VECTOR2D = TransformSteps::Id::Bwd3D::NONE;

            /// Transform step for the third component of a vector
            static const Step VECTOR3D = TransformSteps::Id::Bwd3D::NONE;


            /// Transform step for the first component of a curl
            static const Step CURL1D = TransformSteps::Id::Bwd3D::NONE;

            /// Transform step for the second component of a curl
            static const Step CURL2D = TransformSteps::Id::Bwd3D::NONE;

            /// Transform step for the third component of a curl
            static const Step CURL3D = TransformSteps::Id::Bwd3D::NONE;
         };
      }
   }
}

#endif // TRANSFORMSTEPS_HPP
