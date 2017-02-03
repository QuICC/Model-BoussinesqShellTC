/** 
 * @file SphereTransformSteps.hpp
 * @brief Definition of some useful enums for steps involved in the physical <-> spectral transforms in a whole sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHERETRANSFORMSTEPS_HPP
#define SPHERETRANSFORMSTEPS_HPP

// Configuration includes
//
#include "TypeSelectors/TransformSelector.hpp"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "TransformConfigurators/TransformPath.hpp"

namespace QuICC {

   namespace Transform {

      /**
       * @brief Struct to contain the steps involved in the physical <-> spectral transforms in a whole sphere
       */
      namespace TransformSteps
      {
         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::ProjectorType  Projector1DType;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::ProjectorType  Projector2DType;
         typedef TransformSelector<Dimensions::Transform::TRAND>::Type::ProjectorType  ProjectorNDType;

         typedef TransformSelector<Dimensions::Transform::TRA1D>::Type::IntegratorType  Integrator1DType;
         typedef TransformSelector<Dimensions::Transform::TRA2D>::Type::IntegratorType  Integrator2DType;
         typedef TransformSelector<Dimensions::Transform::TRAND>::Type::IntegratorType  IntegratorNDType;

         /**
          * @brief Generate the list of branches in scalar integration transform
          *
          * @brief isNL Flag to either compute NL projection or pure Physical <-> Spectral
          */
         std::vector<TransformPath>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL);

         /**
          * @brief Generate the list of branches in vector integration transform
          *
          * @brief components Spectral components where to store results: 0: \f$\vec r \nabla\wedge N\f$ 1: \f$\vec r \nabla\wedge\nabla\wedge N\f$
          * @brief isNL Flag to either compute NL projection or pure Physical <-> Spectral
          */
         std::vector<TransformPath>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL);

         /**
          * @brief Generate the list of branches in scalar projection transform
          */
         std::vector<TransformPath>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in scalar gradient transform
          */
         std::vector<TransformPath>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in scalar 2nd order gradient transform
          */
         std::vector<TransformPath>  backwardGradient2(const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req);

         /**
          * @brief Generate the list of branches in vector projection transform
          */
         std::vector<TransformPath>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in vector gradient transform
          */
         std::vector<TransformPath>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in vector curl transform
          */
         std::vector<TransformPath>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req);

         /**
          * @brief Generate the list of branches in vector divergence transform
          */
         std::vector<TransformPath>  backwardDivergence();
      }
   }
}

#endif // SPHERETRANSFORMSTEPS_HPP
