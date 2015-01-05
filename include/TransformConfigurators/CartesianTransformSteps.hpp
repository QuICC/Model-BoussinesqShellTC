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
          * @brief Generate the list of branches in scalar nonlinear term integration transform
          */
         std::vector<IntegratorBranch>  forwardScalarNL();

         /**
          * @brief Generate the list of branches in vector integration transform
          */
         std::vector<IntegratorBranch>  forwardVector();

         /**
          * @brief Generate the list of branches in vector nonlinear term integration transform
          */
         std::vector<IntegratorBranch>  forwardVectorNL();

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
      }
   }
}

#endif // CARTESIANTRANSFORMSTEPS_HPP
