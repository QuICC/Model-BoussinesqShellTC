/** 
 * @file CartesianTransformSteps.cpp
 * @brief Source of the implementation of the cartesian transform steps
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
#include "TransformConfigurators/CartesianTransformSteps.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

namespace TransformSteps {

   std::vector<ProjectorBranch>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::SCALAR));

      return transform;
   }

   std::vector<ProjectorBranch>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::THREE, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::ONE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TWO, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THREE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TWO, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::ONE, FieldType::CURL));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THREE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::ONE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::TWO, FieldType::CURL));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THREE, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::ONE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::CURL));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TWO, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::CURL));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardDivergence()
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::ONE, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::TWO, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::THREE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }

}
}
}
