/** 
 * @file Cartesian2DTransformSteps.hpp
 * @brief Definition of some useful enums for steps involved in a full 2D cartesian transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN2DTRANSFORMSTEPS_HPP
#define CARTESIAN2DTRANSFORMSTEPS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "TransformConfigurators/TransformPath.hpp"

namespace GeoMHDiSCC {

   namespace Transform {

      /**
       * @brief Struct to contain the full field transform steps
       */
      namespace TransformSteps
      {
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

#endif // CARTESIAN2DTRANSFORMSTEPS_HPP