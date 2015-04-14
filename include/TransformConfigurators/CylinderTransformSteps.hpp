/** 
 * @file CylinderTransformSteps.hpp
 * @brief Definition of some useful enums for the steps involved in the physical <-> spectral transforms in a whole cylinder
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CYLINDERTRANSFORMSTEPS_HPP
#define CYLINDERTRANSFORMSTEPS_HPP

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
       * @brief Struct to contain the steps involved in a the physical <-> spectral transforms in a whole cylinder
       */
      namespace TransformSteps
      {
         /**
          * @brief Generate the list of branches in scalar integration transform
          *
          * @brief isNL Flag to either compute NL projection or pure Physical <-> Spectral
          */
         std::vector<IntegratorBranch>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL);

         /**
          * @brief Generate the list of branches in vector integration transform
          *
          * @brief components Spectral components where to store results: 0: \f$\vec r \nabla\wedge N\f$ 1: \f$\vec r \nabla\wedge\nabla\wedge N\f$
          * @brief isNL Flag to either compute NL projection or pure Physical <-> Spectral
          */
         std::vector<IntegratorBranch>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL);

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

#endif // CYLINDERTRANSFORMSTEPS_HPP
